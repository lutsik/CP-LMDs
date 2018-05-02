#
# Very rough code for calling PMDs from window-summarized methylation data using package changepoint
#
#################################### PRELIMINARIES  #################################### 

#### install package changepoint
library(changepoint)
### RnBeads should also be installed, although it is used solely for the methylation color scheme

##### prepare input data first
mmatr  ## a numeric matrix of methylation values (windows x samples). 20kb windows should work O.K.
coords ## a coherent data frame of window genomic coordinates (chr, start, end)

#### set a working directory
ANALYSIS_DIR<-"/tmp"

#################################### START  #################################### 

breakpoints_universe<-list()

ANALYSIS_NAME<-"All_samples_changepoint_pelt_mbic_20kb_final"
dir.create(file.path(ANALYSIS_DIR, ANALYSIS_NAME))
breakpoints_universe[[ANALYSIS_NAME]]<-list()

sample_names<-colnames(mmatr)

for(sample in sample_names){
    
    sample_profile<-mmatr[,samples]
    
    nna.pos<-which(!is.na(sample_profile))
    
    jumps<-which(diff(c(1,nna.pos)) != 1) 
    
    coords.na<-cbind(coords, 
            {v<-rep("", nrow(coords));v[-nna.pos]<-"NA";v},
            {v<-rep("", nrow(coords));v[nna.pos[jumps]]<-"start";v},
            {v<-rep("", nrow(coords));v[nna.pos[jumps-1]]<-"end";v},
            {v<-rep(0, nrow(coords));v[nna.pos[2:length(nna.pos)]]<-diff(nna.pos);v}
    )
    head(coords.na, n=200)
    
    
    intervals.na<-lapply(1:(length(jumps)-1), function(idx) nna.pos[jumps[idx]]:nna.pos[(jumps[idx+1]-1)])
    
    intervals.na.starts<-sapply(1:(length(jumps)-1), function(idx) nna.pos[jumps[idx]])
    intervals.na.ends<-sapply(1:(length(jumps)-1), function(idx) nna.pos[(jumps[idx+1]-1)])
    
    #intervals<-intervals.na
    
    chromosome.starts<-which(diff(as.numeric(coords$V1))>0)
    
    diff_mat_start<-sign(sapply(intervals.na.starts, "-", chromosome.starts))
    diff_mat_end<-sign(sapply(intervals.na.ends, "-", chromosome.starts))
    
    diff_mat_sum<-diff_mat_start+diff_mat_end
    ints_on_chr_border<-apply(diff_mat_sum==0, 1, which)
    
    if(length(ints_on_chr_border)>0){
        for(i in seq_len(length(chromosome.starts))){
            if(length(ints_on_chr_border[[i]])>0){
                intervals.na[[ints_on_chr_border[[i]]]]<-intervals.na.starts[ints_on_chr_border[[i]]]:chromosome.starts[i]
                intervals.na[[length(intervals.na)+1]]<-(chromosome.starts[i]+1):intervals.na.ends[ints_on_chr_border[[i]]]
            }
        }
        intervals.na.starts.fact<-sapply(intervals.na, el, 1)
        
        intervals<-intervals.na[order(intervals.na.starts.fact)]
    }else{
        intervals<-intervals.na
    }
    
    ######## calling
    
    ### changepoint algorithm settings
    
    approach<-c("cpt.mean", "cpt.meanvar")[2]
    
    ### for trying different approaches out
    if(approach=="cpt.meanvar"){
        methods<-c("PELT", "SegNeigh", "BinSeg")[-2]#[c(1,3)]
        penalties<-c( "None", "SIC", "BIC", "MBIC", "AIC", "Hannan-Quinn", "Asymptotic", "Manual", "CROPS")#[4]#[-c(9)]
        default.penalties<-list(0,0,0,0,0,0,0.05,0,c(0,0.1))
        test.stats<-c("Normal","Gamma", "Exponential","Poisson")[-c(2,3,4)]
    }else{
        methods<-c("SegNeigh", "BinSeg")[-1]#[c(1,3)]
        penalties<-c( "None", "SIC", "BIC", "MBIC", "AIC", "Hannan-Quinn", "Asymptotic", "Manual", "CROPS")#[4]#[-c(9)]
        default.penalties<-list(0,0,0,0,0,0,0.05,0,c(0,0.1))
        test.stats<-c("Normal", "CUSUM")
    }
    
    ### for using a selected option 
    methods<-"PELT"
    #methods<-"SegNeigh"
    penalties<-"MBIC"
    #penalties<-"SIC"
    #default.penalties<-list(c(10))
    #default.penalties<-list(c(0,1))
    test.stats<-"Normal"
    
    
    
    ### main loop
    
    for(method.i in 1:length(methods)){
        
        for(penalty.i in 1:length(penalties)){
            if(methods[method.i]=="SegNeigh" && penalties[penalty.i]=="MBIC") next;
            for(ts in 1:length(test.stats)){
                if(test.stats[ts]=="CUSUM" && penalties[penalty.i]%in%c("MBIC", "Asymptotic")) next;
                
                cpt.result<-list()
                for(ii in seq_along(intervals)){
                    if(length(intervals[[ii]])>5){
                        
                        arguments<-list(
                                #data=sample_profile[intervals[[ii]]],
                                data=sample_profile[intervals[[ii]]],
                                #method="BinSeg",
                                method=methods[method.i],
                                test.stat=test.stats[ts],
                                penalty=penalties[penalty.i],
                                Q=length(intervals[[ii]])/2,
                                #Q=length(intervals[[ii]])-1,
                                class=TRUE,
                                param.estimates=TRUE)
                        
                        #arguments$pen.value=default.penalties[[penalty.i]]
                        if(approach=="cpt.meanvar") arguments$shape=1
                        
                        if(penalties[penalty.i]=="Asymptotic") arguments$pen.value=default.penalties[[penalty.i]]
                        if(penalties[penalty.i]=="CROPS") arguments$pen.value=default.penalties[[penalty.i]]
                        if(methods[method.i] %in% c("PELT", "BinSeg")) arguments$minseglen=3
                        
                        cpt.result[[length(cpt.result)+1]]<-do.call(approach, arguments)
                    }else{
                        cpt.result[[length(cpt.result)+1]]<-list()
                        
                    }
                }
                
                #incl<-sapply(cpt.result,length)>0
                #cpt.result<-cpt.result[incl]
                #int.include<-intervals[incl]
                int.include<-intervals
                
                cpt.intervals<-list()
                borders<-list()
                means<-list()
                vars<-list()
                
                for(i in seq_along(cpt.result)){
                    
                    if(length(cpt.result[[i]])>0){
                        cp_list<-cpts(cpt.result[[i]])
                        if(length(cp_list)>0){
                            cpt.intervals[[i]]<-int.include[[i]][cp_list]
                            borders[[i]]<-c(int.include[[i]][1]-1,cpt.intervals[[i]],int.include[[i]][length(int.include[[i]])]+1)
                        }else{
                            cpt.intervals[[i]]<-list()
                            borders[[i]]<-c(int.include[[i]][1]-1,int.include[[i]][length(int.include[[i]])]+1)
                            
                        }
                        means[[i]]<-param.est(cpt.result[[i]])$mean
                        if(approach=="cpt.meanvar") vars[[i]]<-param.est(cpt.result[[i]])$variance
                        
                    }else{
                        borders[[i]]<-c(int.include[[i]][1]-1,int.include[[i]][length(int.include[[i]])]+1)
                        means[[i]]<-mean(sample_profile[int.include[[i]][1]:int.include[[i]][length(int.include[[i]])]])
                        if(approach=="cpt.meanvar") vars[[i]]<-var(sample_profile[int.include[[i]][1]:int.include[[i]][length(int.include[[i]])]])
                    }
                }
                
                ### write out only the called breakpoints
                
                non_zero<-which(sapply(cpt.intervals, length)>0)
                
                cpt.intervals<-unlist(cpt.intervals[non_zero])
                
                cpt.coords<-cbind(coords[cpt.intervals,],cpt.intervals, rep(".", length(cpt.intervals)),  rep(".", length(cpt.intervals)))
                
                write.table(cpt.coords, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("changepoint.coordinates_%s_%s_%s_%s_%s.named.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                cpt.coords[,1]<-gsub("^", "chr", cpt.coords[,1])
                write.table(cpt.coords[,c(1:3)], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("changepoint.coordinates.%s_%s_%s_%s_%s.short.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                breakpoints_universe[[ANALYSIS_NAME]][[sample_names]]<-cpt.intervals
                
                ### write out the complete set of segments
                
                regions<-lapply(seq.int(1,length(borders)), 
                        function(bi) lapply(seq.int(1, length(borders[[bi]])-1), function(bii) (borders[[bi]][bii]+1):(borders[[bi]][bii+1]-1)))
                
                
                regions.starts<-lapply(seq.int(1,length(borders)), 
                        function(bi) lapply(seq.int(1, length(borders[[bi]])-1), function(bii) (borders[[bi]][bii]+1)))
                regions.starts<-unlist(regions.starts)
                
                regions.ends<-lapply(seq.int(1,length(borders)), 
                        function(bi) lapply(seq.int(1, length(borders[[bi]])-1), function(bii) (borders[[bi]][bii+1]-1)))
                regions.ends<-unlist(regions.ends)
                
                seg.means<-unlist(means)
                if(approach=="cpt.meanvar") seg.vars<-unlist(vars) else seg.vars<-rep(0, length(seg.means))
                seg.vars[is.na(seg.vars)]<-0
                
                region_set<-data.frame(Start=regions.starts, End=regions.ends, Mean=seg.means, Var=seg.vars)
                seg.lens<-20*(regions.ends-regions.starts)
                
                filter_low<-seg.means<0.75 & sqrt(seg.vars)<=0.07
                filter_med<-seg.means>=0.5 & sqrt(seg.vars)>=0.07
                filter_high<-seg.means>=0.75 & sqrt(seg.vars)<0.07
                
                region.coord.chroms<-sapply(regions.starts, function(start) coords[start,1])
                region.coord.starts<-sapply(regions.starts, function(start) coords[start,2])
                region.coord.ends<-sapply(regions.ends, function(end) coords[end,3])
                
                
                region.coords<-cbind(gsub("^", "chr", region.coord.chroms), 
                        region.coord.starts, region.coord.ends, rep(".", length(region.coord.chroms)),
                        seg.means, seg.vars)
                
                region.coords[is.na(region.coords[,6]),6]<-0
                
                
                write.table(region.coords, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.coordinates.%s_%s_%s_%s_%s.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                write.table(region.coords[filter_low, ], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.coordinates.%s_%s_%s_%s_%s_low.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                write.table(region.coords[filter_med,], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.coordinates.%s_%s_%s_%s_%s_med.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                write.table(region.coords[filter_high,], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.coordinates.%s_%s_%s_%s_%s_high.bed", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                write.table(region.coords[,c(1:3,5)], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.means.%s_%s_%s_%s_%s.bedGraph", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                write.table(region.coords[,c(1:3,6)], file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, 
                                sprintf("region.vars.%s_%s_%s_%s_%s.bedGraph", paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])), 
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_sd_%s_%s_%s_%s_%s_new.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                plot(seg.means, sqrt(seg.vars), cex=0.1, ylim=c(0,0.4))
                dev.off()
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_sd_smoothed_%s_%s_%s_%s_%s.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                smoothScatter(seg.means, sqrt(seg.vars))
                dev.off()
                
                
                var.bin<-hexbin(seg.means, seg.vars, xbins=75)
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_sd_hexbin_%s_%s_%s_%s_%s.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
                #my_colors=colorRampPalette(c("blue", "yellow", "red"), space = "Lab")
                par(mgp=c(3, 1, 0), cex.axis=1.5)
                plot(var.bin, main="" , colramp=my_colors , trans=log10,
                        #colorcut = seq(0, 1, length=100),
                        legend=F , xlab="segment methylation", ylab="segment variance", clip="off")
                
                dev.off()
                
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_length_smoothed_%s_%s_%s_%s_%s.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                smoothScatter(seg.lens, seg.means)
                dev.off()
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_length_smoothed_zoomed_%s_%s_%s_%s_%s.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                smoothScatter(seg.lens, seg.means, xlim=c(0,100))
                dev.off()
                
                
                len.bin<-hexbin(seg.lens, seg.means, xbins=75)
                
                pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("mean_vs_length_hexbin_%s_%s_%s_%s_%s.pdf",paste(sample_names, collapse="_"), approach, methods[method.i],penalties[penalty.i], test.stats[ts])))
                my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
                #my_colors=colorRampPalette(c("blue", "yellow", "red"), space = "Lab")
                par(mgp=c(3, 1, 0), cex.axis=1.5)
                plot(len.bin, main="" , colramp=my_colors , trans=log10,
                        #colorcut = seq(0, 1, length=100),
                        legend=F , xlab="segment length, kb", ylab="methylation level", clip="off")
                dev.off()
                
            }
        }
    }
}


############################### compare the breakpoints across samples #######################

all.bpts<-sort(unique(Reduce("c", breakpoints_universe[[ANALYSIS_NAME]])))

full_genome<-read.table(file.path(ANALYSIS_DIR, "tracks", "genome.file"))
seqinfo.obj<-Seqinfo(seqnames=as.character(full_genome$V1), seqlengths=full_genome$V2, genome="hg19")
seqinfo.obj<-NULL

granges.all.bpts<-GRanges(coords[all.bpts,1], IRanges(coords[all.bpts,2], coords[all.bpts,3]), "*", seqinfo=seqinfo.obj)

granges.all.reduced<-reduce(granges.all.bpts)


true_bps<-lapply(breakpoints_universe[[ANALYSIS_NAME]], function(bps){
            
            present.bps<-GRanges(coords[bps,1], IRanges(coords[bps,2], coords[bps,3]), "*")
            olaps<-findOverlaps(present.bps, granges.all.reduced)
            unique(subjectHits(olaps))
        })


common.all<-Reduce("intersect", true_bps[names(breakpoints_universe[[ANALYSIS_NAME]])])

pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME,"breakpoints_merged_numbers_pure_more_breakpoints.pdf"), width=5, height=5)
par(mar=c(5.1,4.1,4.1,7.1))
par(oma=c(2,0,0,0))
barplot(sapply(true_bps, length), col=c("salmon"), las=2, ylab="# changepoints")
abline(h=length(common.all), lty=2)
text(14, length(common.all), "common all", cex=0.8, xpd=TRUE, adj=0)
text(14,max(sapply(true_bps, length)), sprintf("total changepoints: %d", length(all.bpts)), xpd=TRUE, cex=0.8)
dev.off()

toBED<-function(granges.obj, rnames=c(rep(".", length(granges.obj))), rscores=c(rep(".", length(granges.obj)))){
    
    df <- data.frame(seqnames=seqnames(granges.obj),
            starts=start(granges.obj)-1,
            ends=end(granges.obj),
            names=rnames,
            scores=rscores,
            strands=strand(granges.obj))
    
    df$starts<-as.integer(df$starts)
    df$ends<-as.integer(df$ends)
    df
}

write.table(toBED(granges.all.reduced), file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"merged_breakpoints.bed"), quote=F, sep="\t", row.names=F, col.names=F)

####### convert breakpoints into methylation domain segments

input.data.coords<-GRanges(coords[,1], IRanges(coords[,2],coords[,3]), strand="*")

granges_full_starts<-GRanges(full_genome$V1, IRanges(1, 2), "*")
granges_full_ends<-GRanges(full_genome$V1, IRanges(full_genome$V2-1,full_genome$V2), "*")

granges.all.reduced.se<-reduce(c(granges_full_starts, granges.all.reduced, granges_full_ends))


olaps<-findOverlaps(input.data.coords,granges.all.reduced.se)

summarized.data.bps<-do.call("rbind", lapply(seq_along(granges.all.reduced.se), function(si) colMeans(mmatr[queryHits(olaps)[subjectHits(olaps)==si],,drop=FALSE])))
nna<-which(rowSums(is.na(summarized.data.bps))==0)

phr<-pheatmap(summarized.data.bps[nna,-10], 
        annotation_col=sample_sheet_wgbs[-10,c(6,8:12)],
        col=RnBeads:::get.methylation.color.panel(),
        scale="none", cluster_true=FALSE,
        cutree_rows=7,
        filename=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"changepoints_heatmap.pdf"))

segments<-gaps(granges.all.reduced.se)

write.table(toBED(segments), file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"merged_segments.bed"), quote=F, sep="\t", row.names=F, col.names=F)


true.seg.lenths<-width(segments)
pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_length_dist.pdf"))
hist(true.seg.lenths, breaks=200,col="lightblue")
dev.off()

##############################################################################

olaps<-findOverlaps(input.data.coords, segments)

summarized.data.segments<-do.call("rbind", lapply(seq_along(segments), function(si) colMeans(mmatr[queryHits(olaps)[subjectHits(olaps)==si],,drop=FALSE], na.rm=TRUE)))
nna<-which(rowSums(is.na(summarized.data.segments))==0)

phr<-pheatmap(summarized.data.segments[nna,-10], 
        annotation_col=sample_sheet_wgbs[,c(6,8:12)],
        col=RnBeads:::get.methylation.color.panel(),
        scale="none", cluster_true=FALSE,
        cutree_rows=6,
        filename=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_heatmap.png"))

domain_clusters<-cutree(phr$tree_row, k=6)
clusters_all<-rep(0, length(segments))
clusters_all[nna]<-as.character(domain_clusters)

##### length distribution by cluster

length.dist<-width(segments)

ldist<-data.frame(Size=length.dist, Log10_size=log10(length.dist), Cluster=clusters_all)
pdf(file.path(ANALYSIS_DIR, ANALYSIS_NAME, "merged_cluster_length_density.pdf"))
ggplot(ldist, aes(x = Log10_size, fill = Cluster)) + geom_density(alpha = 0.5)
dev.off()

###

seg_clusters<-list()
for(mmal in unique(clusters_all)){
    seg_cluster<-segments[clusters_all==mmal]
    seg_clusters[[mmal]]<-seg_cluster
    write.table(toBED(seg_cluster), 
            file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("merged_breakpoints_cluster%s.bed",mmal)), 
            quote=F, sep="\t", row.names=F, col.names=F)
}

### save the first version of the segments
saveRDS(seg_clusters, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_initial.RDS"))

reduced_seg_clusters<-lapply(seg_clusters, function(sc){
            reduce(sc+20000)-20000
        })

### save the reduced version of the segments
saveRDS(reduced_seg_clusters, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_clean.RDS"))


for(rsci in names(reduced_seg_clusters)){
    write.table(toBED(reduced_seg_clusters[[rsci]]), 
            file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("merged_reduced_breakpoints_cluster%s.bed",rsci)), 
            quote=F, sep="\t", row.names=F, col.names=F)
}

################## Refinement of the obtained methylation domains

###### load the cleaned clusters

reduced_seg_clusters<-readRDS(file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_clean.RDS"))
#
##### recluster the merged DMRs
reduced_segments<-Reduce("c", reduced_seg_clusters[as.character(1:6)])

olaps2<-findOverlaps(input.data.coords, reduced_segments)

summarized.data.segments.reduced<-do.call("rbind", lapply(seq_along(reduced_segments), function(si) colMeans(mmatr[queryHits(olaps2)[subjectHits(olaps2)==si],,drop=FALSE], na.rm=TRUE)))
nna<-which(rowSums(is.na(summarized.data.segments.reduced))==0)

##### summarized methylation data of the segments
pdata<-summarized.data.segments.reduced[nna,]
rownames(pdata)<-paste("domain", 1:nrow(pdata))

mcols(reduced_segments)<-pdata

png(file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_heatmap_reduced_new_check.png"), width = 500, height = 500)
phr<-pheatmap:::pheatmap(pdata, 
        annotation_col=sample_sheet_wgbs[,c(6,8:12)],
        col=RnBeads:::get.methylation.color.panel(),
        scale="none", cluster_true=FALSE,
        cutree_rows=6,
        filename=NA)
grid::grid.draw(phr$gtable)
dev.off()

domain_clusters2<-cutree(phr$tree_row, k=6)
clusters_all2<-rep(0, length(reduced_segments))
clusters_all2[nna]<-as.character(domain_clusters2)

#reduced_segments$Cluster<-clusters_all2

##### horizontal heatmap

png(file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_heatmap_reduced_new_horiz.png"), width = 2000, height = 200)
res<-pheatmap:::pheatmap(t(pdata), 
        col=RnBeads:::get.methylation.color.panel(),
        scale="none", cluster_true=FALSE, cutree_cols=6,
        annotation_row=sample_sheet_wgbs[,c(8),drop=FALSE],
        annotation_names_row=FALSE,
        filename=NA
)
grid::grid.draw(res$gtable)
dev.off()

png(file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_heatmap_reduced_new_with_clusters_check.png"))
phr_dummy<-pheatmap:::pheatmap(pdata, 
        annotation_col=sample_sheet_wgbs[,c(6,8:12)],
        annotation_row=data.frame(Cluster=domain_clusters2, row.names=paste("domain", 1:nrow(pdata))), 
        show_rownames=FALSE,
        col=RnBeads:::get.methylation.color.panel(),
        scale="none", cluster_true=FALSE,
        cutree_rows=6,
        filename=NA)
grid::grid.draw(phr_dummy$gtable)
dev.off()

####

seg_clusters_final<-reduced_seg_clusters_refined<-list()

cluster_means<-list()
for(mmal in unique(clusters_all2)){
    cluster_means[[mmal]]<-mean(colMeans(pdata[clusters_all2==mmal,]))
    seg_cluster<-reduced_segments[clusters_all2==mmal]
    reduced_seg_clusters_refined[[mmal]]<-seg_cluster
    seg_clusters_final[[mmal]]<-seg_cluster
}

#### reorder by average methylation
cluster_reorder<-order(unlist(cluster_means)*(sapply(seg_clusters_final,length)>100), decreasing=TRUE)
seg_clusters_final<-seg_clusters_final[cluster_reorder]
reduced_seg_clusters_refined<-reduced_seg_clusters_refined[cluster_reorder]

for(mmal in unique(clusters_all2)){
    write.table(toBED(seg_cluster), 
            file=file.path(ANALYSIS_DIR, ANALYSIS_NAME, sprintf("merged_breakpoints_cluster%s_final.bed",mmal)), 
            quote=F, sep="\t", row.names=F, col.names=F)
}

#### give the segment clusters nice names
names(reduced_seg_clusters_refined)<-as.character(as.roman(seq(length(reduced_seg_clusters_refined))))

#### save the final results
saveRDS(reduced_seg_clusters_refined, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_clean_refined.RDS"))
saveRDS(pdata, file=file.path(ANALYSIS_DIR, ANALYSIS_NAME,"segments_clean_data.RDS"))

#################################### END #################################### 
