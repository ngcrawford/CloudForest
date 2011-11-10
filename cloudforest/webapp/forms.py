from flaskext.wtf import Form, TextField, SelectField, BooleanField, SubmitField, IntegerField
from flaskext.wtf import Required, ValidationError

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