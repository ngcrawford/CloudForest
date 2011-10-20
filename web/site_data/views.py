from site_data import app
from flask import Flask, request, flash, redirect, url_for, render_template
from flaskext.wtf import Form, TextField, Required, BooleanField, SelectField,  RadioField, SubmitField, IntegerField, validators

# app = Flask(__name__)

app.config.update(
    DEBUG=True,
    SECRET_KEY='...'
)

class MyForm(Form):
	aws_id = TextField('AWS ID', validators=[Required()])
	aws_secret_key = TextField('AWS Secret Key', validators=[Required()])
	aws_keypair_name = TextField('AWS Keypair Name', validators=[Required()])
	input_url = TextField('Input Url', validators=[Required()])
	output_url = TextField('Output Url', validators=[Required()])
	job_name = TextField('Job Name', default="CloudForest")
	job_type = SelectField('Job Type', choices=[('genetrees', 'Gene Trees'), ('bootstraps','Bootstraps')])
	mraic = BooleanField('Use MrAIC to infer models',default=True)
	
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






    