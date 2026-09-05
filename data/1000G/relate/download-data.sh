wget https://zenodo.org/records/3234689/files/allele_ages_AFR.zip
unzip allele_ages_AFR.zip
rm allele_ages_AFR.zip

wget https://zenodo.org/records/3234689/files/allele_ages_EUR.zip
unzip allele_ages_EUR.zip
rm allele_ages_EUR.zip

wget https://www.dropbox.com/sh/jmxxp2om4nh3p7u/AADsAtJ5ke8rs4bBHZE6ZDy7a/results/mut
mv mut mut.zip
unzip mut.zip
rm mut.zip
mkdir mut
mv 1000GP* mut/
