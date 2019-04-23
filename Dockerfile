# Build as jupyterhub/singleuser
# Run with the DockerSpawner in JupyterHub

FROM idaes/docker_base:0.1.1

MAINTAINER Project IDAES <scholia@lbl.gov>
USER root

ADD . /home/jovyan/models
RUN chown -R jovyan.users /home/jovyan

USER jovyan

WORKDIR /home/jovyan/models

RUN bash -c 'source activate python2 && conda install -y -c dsdale24 pyqt5'
RUN bash -c 'source activate python2 && pip install -r requirements.txt'
RUN bash -c 'source activate python2 && pyomo --version'
RUN bash -c 'source activate python2 && pip install --upgrade git+https://github.com/pyomo/pyomo.git@ec7cf31c2e7c26f8181e427e73c3f2b8401fef3e'
RUN bash -c 'source activate python2 && pyomo --version'
CMD bash -c 'source activate python2 && nosetests --with-timer'
