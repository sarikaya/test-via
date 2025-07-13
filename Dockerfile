LABEL author="onur.yukselen@umassmed.edu"  description="Docker image containing all requirements for the dolphinnext pipeline"
  
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH


RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git libtbb-dev g++

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# configure image
RUN apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install software-properties-common build-essential
RUN apt install -y dirmngr gnupg apt-transport-https
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt install -y r-base

RUN apt-get update && apt-get install -y gcc unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.0.30.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN aws --version

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda update -n base -c defaults conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext/bin:$PATH

RUN pip install multiqc
RUN pip install wheel
RUN pip install pandas
ENV PATH /usr/local/bin:$PATH
RUN cd /usr/local/bin && wget http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/binaries/reaper-13-100/linux/tally && chmod 777 tally
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
