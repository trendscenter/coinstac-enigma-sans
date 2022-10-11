FROM r-base:4.1.2

# Set the working directory
WORKDIR /computation
RUN apt-get update \
       && apt-get install -y \
       libcurl4-openssl-dev \
       libssl-dev \
       libxml2-dev \
       python3.6 \
       python3-pip \
       vim-tiny \
       && rm -rf /var/lib/apt/lists/*
# Copy requirements.txt into the container
COPY requirements.txt /computation

# Install any needed packages specified in requirements.txt
RUN pip3 install -r requirements.txt

RUN Rscript -e "install.packages(c('rjson','emmeans', 'metafor', 'tidyverse'))"

# Copy the current directory contents into the container
COPY . /computation

CMD ["python3", "entry.py"]
