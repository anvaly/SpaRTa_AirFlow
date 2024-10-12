#!/bin/bash

# ============================================================================================================================== 
# Please do not modify the contents of this script, except for the single section to specify the script(s) below
# ==============================================================================================================================

echo "---- RUNNING run.sh";

# arguments that are provided by the user when triggering the dag and then passed by airflow to this script
download_from_s3=$1;
input_path=$2;
output_path=$3;
updated_factor_sheet=$4;

# update output directory persmissions to allow Anna delete access
update_output_dir_permissions() {
    echo "---- Updating group and permissions of output directory";
    chgrp -R 9500 $output_path;
    chmod -R g+rwx $output_path;
}

# function to exit shell script if previous process failed
exit_if_failed() {
    if [ $? != 0 ]
    then
        echo "ERROR: Previous script failed. Stopping execution of run.sh"
        update_output_dir_permissions;
        exit 1
    fi;
}

# if downloading from S3
if [ $download_from_s3 == 'true' ]
then
    # get credentials to access NGS360 s3 bucket
    secret=$(python get_secret.py);
    export AWS_ACCESS_KEY_ID=$(echo $secret | python -c "import sys, json; print(json.load(sys.stdin)['aws_access_key_id'])");
    export AWS_SECRET_ACCESS_KEY=$(echo $secret | python -c "import sys, json; print(json.load(sys.stdin)['aws_secret_access_key'])");
        
    # download data from NGS360 to stash
    echo "---- RUNNING code/R/a_Wrangle/10X_ngs_s3_to_stash.R";
    Rscript code/R/a_Wrangle/10X_ngs_s3_to_stash.R --meta_file $input_path --output_path $output_path;
    exit_if_failed;

    # reset input_path to point to the factor sheet updated with stash locations
    input_path=$updated_factor_sheet
fi;


# =============================================== MODIFY HERE ==================================================================
# Specify the script(s) to run in the array "scripts_to_run" below 
#   - Note: multiple scripts should be seperated by a single white space, not a comma
#   - Example: scripts_to_run=("path/to/first/script" "path/to/second/script")

scripts_to_run=("code/R/a_Wrangle/10X_to_Seurat_ST.R");
# ==============================================================================================================================

# loops through array, executing each script(s)
for i in "${scripts_to_run[@]}"
do
    echo "---- RUNNING $i";
	Rscript "$i" --meta_file $input_path --output_path $output_path;
    exit_if_failed;
done;

update_output_dir_permissions;