# source of URL list: https://epigenomesportal.ca/api/datahub/download?build=2020-10&assembly=4&format=text

# software to convert from bigBed to bed format
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod a+x bigBedToBed

# files to contain list of relevant datasets
rm -f list-blueprint-h3k4me1.txt
rm -f list-blueprint-h3k4me3.txt

# download and convert Blueprint H3K4me1 / H3K4me3 datasets
mkdir -p bed bigBed

while IFS= read -r url; do

    if [[ $url == "https://epigenomesportal.ca/tracks/Blueprint/hg38/"* ]]; then  # if from Blueprint

        # H3K4me1
        if [[ $url == *"H3K4me1"*".bigBed" ]]; then

            # extract dataset name
            dataset=$(echo $url | awk -F'/' '{print $NF}')

            # download
            wget -O bigBed/$dataset $url

            # make read-only
            chmod 444 bigBed/$dataset

            # extract dataset name without extension
            name=$(echo $dataset | awk -F'.bigBed' '{print $1}')

            # convert
            ./bigBedToBed bigBed/$dataset bed/${name}.bed

            # make read-only
            chmod 444 bed/${name}.bed

            # add to list of relevant datasets
            echo $name >> list-blueprint-h3k4me1.txt
            
            continue  # skip to next dataset
        fi

        # H3K4me3
        if [[ $url == *"H3K4me3"*".bigBed" ]]; then

            # extract dataset name
            dataset=$(echo $url | awk -F'/' '{print $NF}')

            # download
            wget -O bigBed/$dataset $url

            # make read-only
            chmod 444 bigBed/$dataset

            # extract dataset name without extension
            name=$(echo $dataset | awk -F'.bigBed' '{print $1}')

            # convert
            ./bigBedToBed bigBed/$dataset bed/${name}.bed

            # make read-only
            chmod 444 bed/${name}.bed

            # add to list of relevant datasets
            echo $name >> list-blueprint-h3k4me3.txt
        fi
    fi
done < list-urls.txt
