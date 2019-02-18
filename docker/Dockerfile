FROM nfcore/base
LABEL Phil Palmer <phil@lifebit.ai>

# RUN conda config --add channels defaults \
#  && conda config --add channels conda-forge \
#  && conda config --add channels bioconda \
#  && conda install -y salmon=0.8.2 fastqc=0.11.5 multiqc=1.1

# RUN apt-get -y install ttf-dejavu 

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/RNASeq-nf-salmon-dseq2/bin:$PATH

 
