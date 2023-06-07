#!/bin/bash

# Exit on unset variables
set -o nounset
# Exit on unhandled failure in pipes
set -o pipefail
# Have functions inherit ERR traps
set -o errtrace

# Print debug message and terminate script on non-zero return codes
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# [1/2] Most versions of Bash read scripts line by line, causing misbehavior if
# the file changes during runtime. The {} forces Bash to read the entire thing
{
	java -XX:+ExitOnOutOfMemoryError "${@}"

	# [2/2] Prevent Bash from reading past this point once script is done
	exit $?
}
