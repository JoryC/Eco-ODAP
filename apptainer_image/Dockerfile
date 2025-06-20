FROM centos:7.6.1810
LABEL maintainer="Jory Curry <jory.curry@ec.gc.ca>"
LABEL version="1.1"

# Set up Miniconda
ENV MINICONDA_PATH=/opt/miniconda3/bin/
ENV PATH=/opt/miniconda3/bin/:${PATH}
ENV LC_ALL=en_US.utf-8
ENV SHELL=/bin/bash

# Debugging echo statement
RUN echo "Fixing CentOS repository issue by switching to vault.centos.org"

# Install required system dependencies
RUN sed -i 's|^mirrorlist=|#mirrorlist=|g' /etc/yum.repos.d/CentOS-* && \
    sed -i 's|^#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-* && \
    echo "Running yum update" && \
    yum update -y && \
    echo "Installing which, libicu, and curl" && \
    yum install -y which libicu.x86_64 libicu-devel.x86_64 curl && \
    echo "Downloading Miniconda installer" && \
    curl -o Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    echo "Making Miniconda installer executable" && \
    chmod +x Miniconda3-latest-Linux-x86_64.sh && \
    echo "Running Miniconda installer" && \
    ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    echo "Updating PATH and initializing Conda" && \
    export PATH=$PATH:/opt/miniconda3/bin/ && \
    conda init -q bash

# Create a directory for environment files
RUN mkdir -p /tmp/conda_environments

# Copy all environment YAML files
COPY conda_environments/*.yml /tmp/conda_environments/

# Create all Conda environments from YAML files
RUN cd /tmp/conda_environments && \
    for env_file in *.yml; do \
        echo "Creating environment from $env_file" && \
        conda env create -f "$env_file" && \
        echo "Environment from $env_file created successfully"; \
    done && \
    conda clean --all -y

# Remove the YAML files after the environments are created
RUN rm -rf /tmp/conda_environments

# Set up the default shell to use Conda
RUN echo "source /opt/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate" >> ~/.bashrc

# Default command
CMD ["/bin/bash"]