#collect deviant lines

path = '/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/snp/cegs/'
cegs = read.csv(paste0(path,'CEGS.216.lines.NO_DPGP4.GATK.SNP.HETS.FILTERED.11-6-13.tab.txt'), sep = '\t')

#CEGS.216.lines.NO_DPGP4.GATK.SNP.HETS.FILTERED.11-6-13.tab.txt

#gt locus coordinates, X chromosome
gtStart = 2316995
gtEnd = 2334996

cegs[cegs$POS == 2325334, "Raleigh_850"]


#hb locus coordinates, 3R chromosome
hbStart = 4517539
hbEnd = 4535540

#kni locus coordinates, 3L chromosome
kniStart = 20682462
kniEnd = 20700463

#Kr locus coordinates, 2R chromosome
KrStart = 21102137
KrEnd = 21120138

#available positions for giant enhancer
availabilityBoundaries = c(3024, 3279, 6002, 6319, 6965, 8713, 8890, 9196, 9914, 10321, 12015, 12335, 13541, 14291, 15265, 16290, 16712, 17079)
gtBoundaries = availabilityBoundaries + gtStart

#gt regulatory region
cegsX = cegs[cegs[,1] == 'X',]

#find correct indices
gtSNPpositions = cegsX$POS

gtSNPindices = c()
for (i in 1:length(gtSNPpositions)){
    for (j in seq(1,(length(gtBoundaries)-1),2)){
        if (cegsX[i,"POS"] >= gtBoundaries[j] & cegsX[i,"POS"] <= gtBoundaries[j+1]){
            gtSNPindices = c(gtSNPindices, i)
            break
        }
    }
}

#find indices to remove
gtIndicesToRemove = c()
for (i in 1:nrow(cegsX)){
    if ((i %in% gtSNPindices) == FALSE){
        gtIndicesToRemove = c(gtIndicesToRemove,i)
    }
}

#remove SNPs from "closed" regions
cegsX = cegsX[-gtIndicesToRemove,]

#make same levels
for (i in 4:ncol(cegsX)){ levels(cegsX[,i]) = c("./.","A/A","C/C","G/G","T/T") }

#count SNP for gt (X chromosome)
countSNP = c()
for (v in 4:ncol(cegsX)){
    same = 0
    dots = 0
    for (z in 1:nrow(cegsX)){
        same = same + (paste(cegsX[z,3],"/",cegsX[z,3], sep = "") == cegsX[z,v])
        if (cegsX[z,v] == "./."){
            dots = dots + 1
        }
    }
    countSNP = c(countSNP, (nrow(cegsX) - same - dots))
}

#find N lines
N = 40
linesX = c()
for (i in 0:N){
    linesX = c(linesX, colnames(cegsX)[which(countSNP %in% sort(countSNP)[length(countSNP)-i])+3])
}
linesX = unique(linesX)

#find most diversed among those 10
totalMaxDifferent = 0
different = 0
mostDiversedLines = c()
for (i in 1:(length(linesX)-1)){
    
    #save row indices of SNP in current lineX
    firstIndices = c()
    secondIndices = c()
    for (k in 1:nrow(cegsX)){
        # if this is a SNP
        if (paste(cegsX[k,3],"/",cegsX[k,3], sep = "") != cegsX[k,which(colnames(cegsX) == linesX[i])]){
            if (cegsX[k,which(colnames(cegsX) == linesX[i])] != "./."){
                firstIndices = c(firstIndices,k)
            } 
        } else {
            secondIndices = c(secondIndices, k)
        }
    }
    
    #for every other line compare SNPs at different indices
    for (j in (i+1):length(linesX)){
        
        #count different SNPs
        different = 0
        for (p in 1:length(secondIndices)){
            if (paste(cegsX[secondIndices[p],3],"/",cegsX[secondIndices[p],3], sep = "") != cegsX[secondIndices[p],which(colnames(cegsX) == linesX[j])]){
                if (cegsX[secondIndices[p],which(colnames(cegsX) == linesX[j])] != "./."){
                    different = different + 1
                }
            }
        }
        #total amount of SNPs for two lines 
        different = different + length(firstIndices)
        if (different > totalMaxDifferent){
            totalMaxDifferent = different
            mostDiversedLines = c(linesX[i],linesX[j])
        }
    }
}

mostDiversedLines
totalMaxDifferent #69

#collect these SNPs
positions = c()
firstLineSNPs = c()
secondLineSNPs = c()

for (i in 1:nrow(cegsX)){
    positions = c(positions, 18000 - (gtEnd - cegsX$POS[i]))
    
    if (paste0(cegsX[i,mostDiversedLines[1]]) == "./."){
        firstLineSNPs = c(firstLineSNPs, paste0(cegsX[i,"REF"]))
    } else {
        firstLineSNPs = c(firstLineSNPs, paste0(cegsX[i,mostDiversedLines[1]]))
    }
    if (paste0(cegsX[i,mostDiversedLines[2]]) == "./."){
        secondLineSNPs = c(secondLineSNPs, paste0(cegsX[i,"REF"]))
    } else {
        secondLineSNPs = c(secondLineSNPs, paste0(cegsX[i,mostDiversedLines[2]]))
    }
}

#write SNPs in file
write(paste("POS", mostDiversedLines[1], mostDiversedLines[2], sep = "\t"), file = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/gtSNPs.txt", sep = '\n', append = FALSE)
for (i in 1:length(positions)){
    write(paste(positions[i],firstLineSNPs[i], secondLineSNPs[i], sep = "\t"), file = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/gtSNPs.txt", sep = '\n', append = TRUE)
}

#read gt target sites
gtTS = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/gtTScoordinates.txt", header = F, sep = ' ')

#save gt indices of SNPs in target sites
gtSNPindicesTS = c()
for (i in 1:nrow(cegsX)){
    pos = 18002 - (gtEnd - cegsX$POS[i]) - 2
    for (j in 1:nrow(gtTS)){
        if (pos >= gtTS$V3[j] & pos <= gtTS$V4[j]){
            gtSNPindicesTS = c(gtSNPindicesTS,i)
            break
        }
    }
}

#save SNPs from target sites
cegsXts = cegsX[gtSNPindicesTS,]

#count SNP for gt (X chromosome)
countSNPts = c()
for (v in 4:ncol(cegsXts)){
    same = 0
    dots = 0
    for (z in 1:nrow(cegsXts)){
        same = same + (paste(cegsXts[z,3],"/",cegsXts[z,3], sep = "") == cegsXts[z,v])
        if (cegsXts[z,v] == "./."){
            dots = dots + 1
        }
    }
    countSNPts = c(countSNPts, (nrow(cegsXts) - same - dots))
}

#find N lines
N = 10
linesXts = c()
for (i in 0:N){
    linesXts = c(linesXts, colnames(cegsXts)[which(countSNPts %in% sort(countSNPts)[length(countSNPts)-i])+3])
}
linesXts = unique(linesXts)

#find most diversed among those N
totalMaxDifferentTS = 0
different = 0
mostDiversedLinesTS = c()
for (i in 1:(length(linesXts)-1)){
    
    #save row indices of SNP in current lineX
    firstIndices = c()
    secondIndices = c()
    for (k in 1:nrow(cegsXts)){
        # if this is a SNP
        if (paste(cegsXts[k,3],"/",cegsXts[k,3], sep = "") != cegsXts[k,which(colnames(cegsXts) == linesXts[i])]){
            if (cegsXts[k,which(colnames(cegsXts) == linesXts[i])] != "./."){
                firstIndices = c(firstIndices,k)
            } else {
                secondIndices = c(secondIndices, k)
            }
        } else {
            secondIndices = c(secondIndices, k)
        }
    }
    
    #for every other line compare SNPs at different indices
    for (j in (i+1):length(linesXts)){
        
        #count different SNPs
        different = 0
        for (p in 1:length(secondIndices)){
            if (paste(cegsXts[secondIndices[p],3],"/",cegsXts[secondIndices[p],3], sep = "") != cegsXts[secondIndices[p],which(colnames(cegsXts) == linesXts[j])]){
                if (cegsXts[secondIndices[p],which(colnames(cegsXts) == linesXts[j])] != "./."){
                    different = different + 1
                }
            }
        }
        #total amount of SNPs for two lines 
        different = different + length(firstIndices)
        if (different > totalMaxDifferentTS){
            totalMaxDifferentTS = different
            mostDiversedLinesTS = c(linesXts[i],linesXts[j])
        }
    }
}

mostDiversedLinesTS
totalMaxDifferentTS #19

#collect these SNPs
positions = c()
firstLineSNPs = c()
secondLineSNPs = c()

for (i in 1:nrow(cegsX)){
    positions = c(positions, 18000 - (gtEnd - cegsX$POS[i]))
    
    if (paste0(cegsX[i,mostDiversedLinesTS[1]]) == "./."){
        firstLineSNPs = c(firstLineSNPs, paste0(cegsX[i,"REF"]))
    } else {
        firstLineSNPs = c(firstLineSNPs, paste0(cegsX[i,mostDiversedLinesTS[1]]))
    }
    if (paste0(cegsX[i,mostDiversedLinesTS[2]]) == "./."){
        secondLineSNPs = c(secondLineSNPs, paste0(cegsX[i,"REF"]))
    } else {
        secondLineSNPs = c(secondLineSNPs, paste0(cegsX[i,mostDiversedLinesTS[2]]))
    }
}

#write SNPs in TS to file
write(paste("isTS", "POS", mostDiversedLinesTS[1], mostDiversedLinesTS[2], sep = "\t"), file = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/gtSNPts.txt", sep = '\n', append = FALSE)
for (i in 1:length(positions)){
    if (positions[i] %in% (18000 - (gtEnd - cegsXts$POS))){
        write(paste(1,positions[i],firstLineSNPs[i], secondLineSNPs[i], sep = "\t"), file = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/gtSNPts.txt", sep = '\n', append = TRUE)
    } else {
        write(paste(0,positions[i],firstLineSNPs[i], secondLineSNPs[i], sep = "\t"), file = "/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/biodata/Drosophila/genes&chromatin/gtSNPts.txt", sep = '\n', append = TRUE)
        
    }
}

#countSNPts[which(colnames(cegsXts) == mostDiversedLinesTS[2])-3]

#___________________________________________________________________________

#open
openSNPs = 0

#save row indices of SNP
firstIndices = c()
secondIndices = c()
for (k in 1:nrow(cegsX)){
    # if this is a SNP
    if (paste(cegsX[k,3],"/",cegsX[k,3], sep = "") != cegsX[k,mostDiversedLinesTS[1]]){
        if (cegsX[k,mostDiversedLinesTS[1]] != "./."){
            firstIndices = c(firstIndices,k)
        } 
    } else {
        secondIndices = c(secondIndices, k)
    }
}
#count different SNPs
for (p in 1:length(secondIndices)){
    if (paste(cegsX[secondIndices[p],3],"/",cegsX[secondIndices[p],3], sep = "") != cegsX[secondIndices[p],mostDiversedLinesTS[2]]){
        if (cegsX[secondIndices[p],mostDiversedLinesTS[2]] != "./."){
            openSNPs = openSNPs + 1
        }
    }
}
#total amount of SNPs for two lines 
openSNPs = openSNPs + length(firstIndices)
openSNPs

#ts
tsSNPs = 0

#save row indices of SNP in current
firstIndices = c()
secondIndices = c()
for (k in 1:nrow(cegsXts)){
    # if this is a SNP
    if (paste(cegsXts[k,3],"/",cegsXts[k,3], sep = "") != cegsXts[k,mostDiversedLines[1]]){
        if (cegsXts[k,mostDiversedLines[1]] != "./."){
            firstIndices = c(firstIndices,k)
        } else {
            secondIndices = c(secondIndices, k)
        }
    } else {
        secondIndices = c(secondIndices, k)
    }
}

#for every other line compare SNPs at different indices
#count different SNPs
for (p in 1:length(secondIndices)){
    if (paste(cegsXts[secondIndices[p],3],"/",cegsXts[secondIndices[p],3], sep = "") != cegsXts[secondIndices[p],mostDiversedLines[2]]){
        if (cegsXts[secondIndices[p],mostDiversedLines[2]] != "./."){
            tsSNPs = tsSNPs + 1
        }
    }
}
#total amount of SNPs for two lines 
tsSNPs = tsSNPs + length(firstIndices)
tsSNPs

#analysis
countSNP[which(colnames(cegsX) == mostDiversedLines[1])-3]
countSNP[which(colnames(cegsX) == mostDiversedLines[2])-3]
countSNPts[which(colnames(cegsX) == mostDiversedLines[1])-3]
countSNPts[which(colnames(cegsX) == mostDiversedLines[2])-3]

countSNP[which(colnames(cegsX) == mostDiversedLinesTS[1])-3]
countSNP[which(colnames(cegsX) == mostDiversedLinesTS[2])-3]
countSNPts[which(colnames(cegsX) == mostDiversedLinesTS[1])-3]
countSNPts[which(colnames(cegsX) == mostDiversedLinesTS[2])-3]


