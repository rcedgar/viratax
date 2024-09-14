#!/bin/bash

sudo docker run --rm --mount type=bind,source=.,target=/docker_results viratax
