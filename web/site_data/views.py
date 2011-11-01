from site_data import app
# from secret_data import fake_form_data
from flask import Flask, request, jsonify, flash, redirect, url_for, render_template
from flaskext.wtf import Form, TextField, Required, BooleanField, SelectField,  RadioField, SubmitField, IntegerField, validators, ValidationError

from CloudForest.cloudforest import cloudforest
from multiprocessing import Process


# app = Flask(__name__)

app.config.update(
    DEBUG=True,
    SECRET_KEY='...'
)


class MyForm(Form):
    """ Create cloudforest form"""
    
    def validate_S3_url(form, field):
        """Checks that field is a valid S3 url.
            Needs improvement. :) """
        if field.data == "s3://":
            raise ValidationError(u'This field must be an S3 url.')

    # Setup fields
    aws_id = TextField('AWS ID', validators=[Required()])
    aws_secret_key = TextField('AWS Secret Key', validators=[Required()])
    aws_region = SelectField('AWS Region', choices=[('us_east_1', 'US East'), ('us_west_1', 'US West')])
    input_url = TextField('Input Url', validators=[Required(), validate_S3_url])
    output_url = TextField('Output Url', validators=[Required(), validate_S3_url])
    job_name = TextField('Job Name', default="CloudForest", validators=[Required()])
    job_type = SelectField('Job Type', choices=[('genetrees', 'Gene Trees'), ('bootstraps','Bootstraps')])
    mraic = BooleanField('Use MrAIC to infer models',default=False)
    map_tasks = IntegerField('Map Tasks', default=19)
    reduce_tasks = IntegerField('Reduce Tasks', default=1)
    bootstraps = IntegerField('Bootstrap Replicates', default=0)
    ec2_instance_type = SelectField('EC2 Instance Type:', choices=[("m1.small","m1.small")])
    
def run_cloudforest(form_options):
    """Take form options and inject them into cloudforest"""
    cf_options = ['-r', 'emr', '--no-conf', \
                '--num-ec2-instances', str(2), \
                '--jobconf', ("aws_access_key_id=%s" % (form_options['aws_id'])), \
                '--jobconf', ("aws_secret_access_key=%s" % (form_options['aws_secret_key'].upper())), \
                '--jobconf', ("aws_region=%s" % (form_options['aws_region'])), \
                '--jobconf', ('mapred.map.tasks=%s' % (form_options['map_tasks'])), \
                '--jobconf', ('mapred.reduce.tasks=%s' % (form_options['reduce_tasks'])), \
                '--jobconf', 'mapred.reduce.tasks.speculative.execution=True', \
                '--archive', 's3://bioaws/mapreduce/exes/aws.phylo.tar.gz#bin',]
    
    if form_options['mraic'] != None:
        cf_options.append("--mraic")
    
    if form_options['job_type'] == 'genetrees':
        cf_options.append("--gene-trees")
    
    if form_options['bootstraps'] != 0:
        cf_options.append("--bootstraps")
        cf_options.append(form_options['bootstraps'])
        
    cf_options.append(form_options['input_url'])
    
    mr_job = cloudforest.ProcessPhyloData()
    mr_job.load_options(cf_options)
    print 'submit executed'
    
    # with mr_job.make_runner() as runner:
    #     runner.run()
    #     for line in runner.stream_output():
    #         key, value = mr_job.parse_output_line(line)
    #         print key, value


@app.route("/", methods=("GET", "POST"))
def submit():
    """Submits form data to cloudforest"""
    form = MyForm(request.form, csrf_enabled=False)
    if form.validate_on_submit(): 
        # I think what we want to do here is spawn off the cloudforest into its 
        # own subprocess. It sort of works for now.
        
        # ===========================
        #  COMMENT OUT LINES 88-90 UNLESS
        #  YOU WANT TO INSTALL A BUNCH
        #  OF PYTHON SCIENTIFIC DEPENDANCIES
        
        p = Process(target=run_cloudforest(form.data))
        p.start()
        p.join()
        
        pass # don't comment out the pass or the if statement becomes invalid
        
    return render_template("index.html",form=form)






    