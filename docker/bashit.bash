#!/bin/bash

sudo docker run --rm -it --mount type=bind,source=.,target=/docker_results viratax bash
