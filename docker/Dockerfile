FROM ubuntu:22.04
RUN apt-get update

# Avoid interactive prompts for geographic area / timezone
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get install -y wget
RUN apt-get install -y python3
RUN apt-get install -y git
RUN apt-get install -y emboss

WORKDIR /root

RUN wget https://serratus-public.s3.amazonaws.com/rce/viratax/viratax-repo-2024-09-14.tar.gz
RUN tar -zxvf viratax-repo-2024-09-14.tar.gz

WORKDIR /root/scripts

CMD [ "/bin/bash", "-c", "./RUNME.bash" ]
