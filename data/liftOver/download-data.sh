# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz

# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
