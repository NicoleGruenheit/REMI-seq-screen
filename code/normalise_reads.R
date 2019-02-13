normalise_reads <- function(df,comb_choice,cutoff) {
  
  # get rid of tags that map to multiple positions
  df = filter(df,!V28 == "multiple")
  df = filter(df,(V10 == "GATC" & (V3 == "G1" | V3 == "G2" | V3 == "G3")) | (V10 == "CATG" & (V3 == "C4" | V3 == "C7" | V3 == "C8")))
  # add info about the possible disrupted gene
  df = mutate(df, dis = ifelse(V24 == "+" & V25 == "+",ifelse(V27 < 500,"intergenic_down","intergenic"),ifelse(V24 == "-" & V25 == "+",ifelse(V27 < 500 & V26 < 500,"intergenic_both","intergenic"),ifelse(V24 == "-" & V25 == "-",ifelse(V26 < 500,"intergenic_up","intergenic"),"intergenic"))))
  df = mutate(df, dis = ifelse(! is.na(V16),"intragenic",dis))
  df = mutate(df, dis_gene = ifelse(dis == "intragenic",V16,ifelse(dis == "intergenic","none",ifelse(dis == "intergenic_down",V23,ifelse(dis=="intergenic_up",V22,ifelse(dis == "intergenic_both",paste(V22,V23,sep="-"),NA))))))
  if(comb_choice==1){
    all = group_by(df,V6,V7,V30)
    all = select(all,-V3)                                   # get rid of the index
    all = mutate(all,tag = sum(V5))                         # add up counts for both sides
    all = unique(select(all,V6,V7,V30,dis,dis_gene,tag))    # loose half of the data
    all = group_by(all,V30)                     # normalise by total number of reads per sample
    all = mutate(all,total=sum(tag),norm = 1000000 * tag/total)
    # loose a lot of the columns because we don't need them for the rest (can add back later)
    all = select(all,V6,V7,V30,dis,dis_gene,norm)
    all = unique(all)
  }else if(comb_choice==2){
    all= group_by(df,V3,V6,V7,V30)
    all = mutate(all,tag = sum(V5))
    all = group_by(all,V30)
    all = mutate(all,total=sum(tag),norm = 1000000 * tag/total)
    # count number of reads per position (grouping by index (V3), name, and sample)
    all = group_by(all,V6,V7,V30)
    all = mutate(all,total_index = 1000000 *sum(V5)/total)
    all = mutate(all,filter = ifelse(total_index < 5 * norm/100,1,0))
    all = filter(all,filter==0)
    all = select(all,V3,V6,V7,V30,dis,dis_gene,norm)
    all = unique(all)
  }else if(comb_choice==3){
    
    df  = filter(df,!is.na(V16))
    all = group_by(df,dis_gene,V30)
    all = select(all,-V3,-V6,-V7)                                   # get rid of the index
    all = mutate(all,tag = sum(V5))
    all = group_by(all,V30)                                         # normalise by total number of reads per sample
    all = mutate(all,total=sum(tag),norm = 1000000 * tag/total)
    # loose a lot of the columns because we don't need them for the rest (can add back later)
    all = select(all,V30,dis,dis_gene,norm)
    all = unique(all)
    print(head(all))
  }else if(comb_choice==4){
    df  = filter(df,!(dis == "intergenic" | dis == "intergenic_both"))
    all = group_by(df,dis_gene,V30)
    all = select(all,-V3,-V6,-V7)                                   # get rid of the index
    all = mutate(all,tag = sum(V5))
    all = group_by(all,V30)                                         # normalise by total number of reads per sample
    all = mutate(all,total=sum(tag),norm = 1000000 * tag/total)
    # loose a lot of the columns because we don't need them for the rest (can add back later)
    all = select(all,V30,dis_gene,norm)
    all = unique(all)
  }else if(comb_choice==5){
    all = group_by(df,V6,V7,V30)
    all = select(all,-V3)                                   # get rid of the index
    all = mutate(all,norm = sum(V5))                         # add up counts for both sides
    all = unique(select(all,V6,V7,V30,dis,dis_gene,norm))    # loose half of the data
    # loose a lot of the columns because we don't need them for the rest (can add back later)
    all = select(all,V6,V7,V30,dis,dis_gene,norm)
    all = unique(all)
  }
  if(!(comb_choice == 3 | comb_choice == 4)){
    all = group_by(all,V6,V7)
    aa = summarise(all,counts = sum(norm > cutoff, na.rm = TRUE))
    print(dim(aa))
    aa = filter(aa,counts > 0)
    aa = mutate(aa,name = paste(V6,V7,sep="_"))
    print(dim(aa))
    all = mutate(all,name = paste(V6,V7,sep="_"))
  }else {
    all = group_by(all,dis_gene)
    aa = summarise(all,counts = sum(norm > cutoff, na.rm = TRUE))
    print(dim(aa))
    aa = filter(aa,counts > 0)
    aa = mutate(aa,name = dis_gene)
    print(dim(aa))
    all = mutate(all,name = dis_gene)
  }
  
  all = filter(all,name %in% aa$name)
  all = select(all,-name)
  all = spread(all,V30,norm)                  # convert into wide format
  all = as.data.frame(all)
  all[is.na(all)] = 1
  all[all < 1] = 1
  all = ungroup(all)
  return(all)
}