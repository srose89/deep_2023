FROM python:3.8

MAINTAINER roses3@mskcc.org

USER root

WORKDIR /app

ADD . /app

ARG HDF5_DIR=/opt/homebrew/opt/hdf5
ARG BLOSC_DIR=/opt/homebrew/opt/c-blosc

RUN pip install --trusted-host pypi.python.org -r requirements.txt

EXPOSE 8050

CMD ["python", "application.py"]
