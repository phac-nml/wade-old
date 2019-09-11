#!/bin/bash 

while true; do
	sleep 10m
	if [ `netstat -t | grep -v CLOSE_WAIT | grep ':80' | wc -l` -lt 3 ]
    then
		pkill shiny-server
	fi
done
