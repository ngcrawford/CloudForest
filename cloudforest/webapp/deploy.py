import cloudforest
from forms import MyForm
from multiprocessing import Process
from flask import Flask, request, render_template

application = Flask(__name__)
application.config.from_object(__name__)

application.config.update(
    DEBUG=True,
    SECRET_KEY='...'
)

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
        cf_options.applicationend("--mraic")
    
    if form_options['job_type'] == 'genetrees':
        cf_options.applicationend("--gene-trees")
    
    if form_options['bootstraps'] != 0:
        cf_options.applicationend("--bootstraps")
        cf_options.applicationend(form_options['bootstraps'])
        
    cf_options.applicationend(form_options['input_url'])
    
    mr_job = cloudforest.ProcessPhyloData()
    mr_job.load_options(cf_options)
    print 'submit executed'
    
    # with mr_job.make_runner() as runner:
    #     runner.run()
    #     for line in runner.stream_output():
    #         key, value = mr_job.parse_output_line(line)
    #         print key, value

@application.route("/", methods=("GET", "POST"))
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

if __name__ == '__main__':
    application.run()
