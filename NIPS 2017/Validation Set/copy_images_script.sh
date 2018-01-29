#!/bin/bash
root_folder=/mnt/disks/DATA/raw_resized/Data/CLS-LOC/train/
folders="$(gcloud compute ssh junjiajia_long@create-adv --command="ls $root_folder")"
for folder in $folders
do
	mkdir -p images/$folder
	sub_folder=$root_folder$folder/
	echo copying from $sub_folder
	files="$(gcloud compute ssh junjiajia_long@create-adv --command="ls $sub_folder | sed -n -e '11,20p'")"
	for file in $files
	do
		source=$sub_folder$file
		dest=images/$folder/$file
		echo copying from $source to $dest
		gcloud compute scp junjiajia_long@create-adv:$source $dest
	done
	echo
done