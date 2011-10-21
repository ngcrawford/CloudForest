from site_data import app
from flask import Flask, request, flash, redirect, url_for, render_template
from flaskext.wtf import Form, TextField, Required, BooleanField, SelectField,  RadioField, SubmitField, IntegerField, validators, ValidationError

# app = Flask(__name__)

app.config.update(
    DEBUG=True,
    SECRET_KEY='...'
)

class MyForm(Form):

	def validate_S3_url(form, field):
		"""Checks that field is a valid S3 url.
			Needs improvement. :)
		"""
		if field.data == "s3n://":
			raise ValidationError(u'This field must be an S3 url.')

	aws_id = TextField('AWS ID', validators=[Required()])
	aws_secret_key = TextField('AWS Secret Key', validators=[Required()])
	aws_keypair_name = TextField('AWS Keypair Name', validators=[Required()])
	input_url = TextField('Input Url', validators=[Required(), validate_S3_url])
	output_url = TextField('Output Url', validators=[Required(), validate_S3_url])
	job_name = TextField('Job Name', default="CloudForest", validators=[Required()])
	job_type = SelectField('Job Type', choices=[('genetrees', 'Gene Trees'), ('bootstraps','Bootstraps')])
	mraic = BooleanField('Use MrAIC to infer models',default=True)
	map_tasks = IntegerField('Map Tasks', default=19)
	reduce_tasks = IntegerField('Reduce Tasks', default=1)
	bootstraps = IntegerField('Bootstrap Replicates', default=0)
	

	
	
	
# if os.sys.platform == 'darwin':
#   result_url = "/query_args"
# else:
#   result_url = "/ngcrawford/query_args"


@app.route("/", methods=("GET", "POST"))
def submit():
	form = MyForm(request.form, csrf_enabled=False)
	if form.validate_on_submit():
	    print form.data

	return render_template("index.html", form=form)






    