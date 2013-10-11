#!/usr/bin/ruby
#
# Short ruby script to watch over the progress of a running
# Nomadm/TLEED calculation
# This script runs the tleed_progress script then
# goes to sleep for 5 minutes. You can adjust how many times you want to do
# this by changing the value of the max count, i.e. 50
# The tlee_progress script checks the tmp.log file in this directory and
# outputs a short summary of the Nomad results to date
#
count = 0
while count < 50
  shell_script_output = `./tleed_progress >> pertsummary` 
#  shell_script_output = `./shell-script.sh`
  sleep 300
  count=count+1
  puts count
end
