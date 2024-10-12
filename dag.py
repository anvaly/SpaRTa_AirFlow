from datetime import datetime, timedelta
import os
from datetime import datetime
from urllib.parse import quote_plus

from airflow import DAG
from airflow.operators.python import PythonOperator, BranchPythonOperator
from airflow.operators.email import EmailOperator
from airflow.operators.bash import BashOperator
from plugins.ecs_ec2_operator import ECSEC2Operator


app_name = 'sparta_automation'
cluster = app_name + 'Cluster'
container_name = app_name + 'Container'
task_definition = app_name + 'TaskDef'
task_execution_role = 'arn:aws:iam::####:role/XXXX'
launch_type = 'EC2'
log_group = app_name + 'LogGroup'

default_args = {
    'owner': 'Daniel Carrera',
    'email': ['XXXX'],
    'email_on_failure': True,
    'email_on_retry': False,
    'start_date': datetime(2023, 3, 30),
    'retries': 0,
    'retry_delay': timedelta(minutes=5),
    'params': {
        'INSTRUCTIONS':'email: must be a BMS email address. download_from_s3: must be either true or false. input_path: path to a factor sheet that contains s3 locations if downloading from s3, or stash paths otherwise. output_folder_name: name of the folder to output to, this name will be concatenated with a timestamp, if no name is desired, then leave blank as is', 
        'email': 'REQUIRED', 
        'download_from_s3': 'REQUIRED',
        'input_path': 'REQUIRED', 
        'output_folder_name': ''
        }
}

network_configuration = {
        'awsvpcConfiguration': {
            'subnets': ['XXXX', 'XXXX'],
            'securityGroups': ['XXXX', 'XXXX']
        }
    }

# needed to mount /XXXX to ec2 instance
user_data = f"""#!/bin/bash
echo ECS_CLUSTER={cluster} >> /etc/ecs/ecs.config
echo ECS_ENABLE_TASK_IAM_ROLE=true >> /etc/ecs/ecs.config
echo ECS_ENABLE_CONTAINER_METADATA=true >> /etc/ecs/ecs.config
echo 'user_allow_other' >> /etc/fuse.conf
sudo amazon-linux-extras install epel -y
sudo yum update -y && sudo yum install -y sshfs jq unzip
sudo curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
sudo unzip awscliv2.zip
sudo ./aws/install
sudo mkdir /XXXX
sudo chmod 777 /XXXX
export AWS_DEFAULT_REGION=us-east-1
/usr/local/bin/aws secretsmanager get-secret-value --secret-id XXXX | jq --raw-output '.SecretString' | jq -r .password | sshfs XXXX@XXXX:/XXXX /XXXX -o password_stdin -o StrictHostKeyChecking=no -o allow_other -o idmap=user
"""


def _get_user_input(ti, **kwargs):
    """
    User provided arguments:
    - email: must be a BMS email address.(REQUIRED)
    - download_from_s3: must be either true or false. 
    - input_path: path to a factor sheet that contains s3 locations if downloading from s3, or stash paths otherwise. (REQUIRED)
    - output_folder_name: name of the folder to output to, this name will be concatenated with a timestamp. (OPTIONAL)
    
    Note: input_path is not checked for existence as this requires the mount as service account user for file read permissions. This mount occurs in the ECS cluster where Anna's R code does these checks. 
    """
    # set email recipients
    try:
        email_recipients = kwargs["dag_run"].conf["email"]
        assert(email_recipients[-7:] == 'bms.com')
    except:
        # send anna email if email not provided
        ti.xcom_push(key='email_recipients', value='XXXX')
        ti.xcom_push(key='error_message', value="ERROR: The Sparta pipeline was triggered without an email or with a non-BMS email. Please provide a BMS email and retrigger.")
        return 'send_user_input_email'
    else:
        ti.xcom_push(key='email_recipients', value=email_recipients)

    # get required arguments
    try:    
        download_from_s3 = kwargs["dag_run"].conf["download_from_s3"]
        assert(download_from_s3.lower() in ['true', 'false'])
        
        input_path = kwargs["dag_run"].conf["input_path"]
        
        output_folder_name = kwargs["dag_run"].conf["output_folder_name"]     
    except:
       ti.xcom_push(key='error_message', value="ERROR: One or more of the required arguments were not provided when triggering the DAG or are not valid values. Please refer to the instructions listed when triggering the DAG.")
       return 'send_user_input_email'
    
    # construct output_path and check if exists - mount as resvc not required as airflow instance has read permissions to this location
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    if output_folder_name:
        output_path = f"/XXXX/{output_folder_name}_{timestamp}"
    else:
        output_path = f"/XXXX/{timestamp}"  
        
    if os.path.exists(output_path):
        ti.xcom_push(key='error_message', value=f"ERROR: The following output path already exists: {output_path}")
        return 'send_user_input_email'

    # construct path to updated factor sheet
    updated_factor_sheet_path = os.path.join(output_path, '1_SpaceRanger_outputs', 'meta_file.txt')

    # get dag run's execution date, which will be used to construct the url to it's log on the web browser
    execution_date = kwargs["ts"]
    execution_date_html = quote_plus(execution_date) # get in required encoding
    
    # # get required instance_type
    # try:
    #     instance_type = kwargs["dag_run"].conf["instance_type"]
    # except:
    #    ti.xcom_push(key='error_message', value="ERROR: The Sparta pipeline was triggered without a specified instance_type. Please retrigger and provide an instance_type.")
    #    return 'send_user_input_email'

    # push xcoms
    ti.xcom_push(key='download_from_s3', value=download_from_s3.lower())
    ti.xcom_push(key='input_path', value=input_path)
    ti.xcom_push(key='updated_factor_sheet_path', value=updated_factor_sheet_path)
    ti.xcom_push(key='timestamp', value=timestamp)
    ti.xcom_push(key='output_path', value=output_path)
    ti.xcom_push(key='execution_date_html', value=execution_date_html) 
    # ti.xcom_push(key='instance_type', value=instance_type) 
    return 'run_etl'


with DAG(
    app_name + 'ETL',
    default_args=default_args,
    schedule_interval=None,
    catchup=False
) as dag:

  get_user_input = BranchPythonOperator(
    task_id = 'get_user_input',
    python_callable = _get_user_input,
    provide_context = True
  )

  send_user_input_email = EmailOperator(
    task_id='send_user_input_email', 
    to= ["{{ task_instance.xcom_pull(key='email_recipients') }}"], 
    subject="""Sparta Automated Pipeline: Error""",  
    html_content="""
      <html>
      <p>Hello,</p>
      <p>
      {{ task_instance.xcom_pull(key='error_message') }}
      </p>
      <br>
      <br>
      If you would like to be removed from this emailing list, please reach out to: XXXX
      <br>
      <p>Thank you,</p>
      </html>
    """
  )
  
  run_etl = ECSEC2Operator(
    task_id='run_etl',
    dag=dag,
    task_definition=task_definition,
    instance_type="r5a.4xlarge",
    cluster=cluster,
    network_configuration=network_configuration,
    container_name=container_name,
    container_command=["/bin/bash", "run.sh", "{{ task_instance.xcom_pull(key='download_from_s3') }}", "{{ task_instance.xcom_pull(key='input_path') }}", "{{ task_instance.xcom_pull(key='output_path') }}", "{{ task_instance.xcom_pull(key='updated_factor_sheet_path') }}"],
    ec2_user_data=user_data,
    key_pair_name='XXXX',
    awslogs_group=log_group,
    awslogs_stream_prefix=app_name,
    ami_id='XXXX',
    sleep_time=240,
    )
  
  success_email = EmailOperator(
    task_id='success_email', 
    to= ["{{ task_instance.xcom_pull(key='email_recipients') }}"], 
    subject="""Sparta Automated Pipeline Completed""",  
    html_content="""
      <html>
      <p>Hello,</p>
      <p>
      The following Sparta automated pipeline run was just completed:
      <a href="XXXX">Click here to view logs!</a>
      <br>
      <br>
        Timestamp: {{ task_instance.xcom_pull(key='timestamp') }}
      <br>
      <br>
        Input factor sheet: {{ task_instance.xcom_pull(key='input_path') }}
      <br>
      <br>
        Data downloaded from NGS360 bucket: {{ task_instance.xcom_pull(key='download_from_s3') }}
      <br>
      <br>
        Output location: {{ task_instance.xcom_pull(key='output_path') }}
      </p>
      <br>
      <br>
      If you would like to be removed from this emailing list, please reach out to: XXXX
      <br>
      <p>Thank you,</p>
      </html>
    """
  )

  fail_email = EmailOperator(
    task_id='fail_email', 
    to= ["{{ task_instance.xcom_pull(key='email_recipients') }}"], 
    subject="""Sparta Automated Pipeline Failed""",  
    html_content="""
      <html>
      <p>Hello,</p>
      <p>
      The following Sparta automated pipeline just failed. Please see logs to view error(s). <br>
      <a href="XXXX">Click here to view logs!</a>
      <br>
      <br>
        Timestamp: {{ task_instance.xcom_pull(key='timestamp') }}
      <br>
      <br>
        Input factor sheet: {{ task_instance.xcom_pull(key='input_path') }}
      <br>
      <br>
        Data downloaded from NGS360 bucket: {{ task_instance.xcom_pull(key='download_from_s3') }}
      <br>
      <br>
        Output location: {{ task_instance.xcom_pull(key='output_path') }}
      </p>
      <br>
      <br>
      If you would like to be removed from this emailing list, please reach out to: XXXX
      <br>
      <p>Thank you,</p>
      </html>
    """,
    trigger_rule = 'one_failed'
  )

  
  get_user_input >> send_user_input_email
  get_user_input >> run_etl >> success_email
  run_etl >> fail_email