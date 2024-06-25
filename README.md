# endoderm-app
Dash app for the mouse endoderm data explorer

### Elastic Beanstalk

Branch `eb-app` is used for Elastic Beanstalk. 


#### Config files
The data used is stored on an EFS volume. Creation and mounting of EFS volumes is performed using the AWS provided config files in `.ebextensions`. 

Python apps also require an additional config file: `.ebextensions/wsgi_custom.config`


#### Deployment
* Install AWS elastic beanstalk CLI (EB CLI): https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/eb-cli3.html
* `eb init` to initialize an application
* `eb create` to create an environment and run everything.
