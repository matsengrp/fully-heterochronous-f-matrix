default:

install:
	pip install -e '.[dev]'

test:
	pytest tests

format:
	docformatter --in-place --black --recursive farmrats tests
	black farmrats tests

checkformat:
	docformatter --check --black --recursive farmrats tests
	black --check farmrats tests

checktodo:
	grep -rq --include=\*.{py} "TODO" . && echo "TODOs found" && exit 1 || echo "No TODOs found" && exit 0

lint:
	flake8 . --max-complexity=30 --ignore=E731,W503,E402,F541,E501,E203,E266 --statistics --exclude=_ignore
# Ignored errors complain about...
# E731: lambda functions
# W503: order of line breaks and operators in multiline commands
# E402: imports not all at top of file
# F541: unnecessary f-strings 
# E501: long lines
# E203: whitespace before comma, colon, or semicolon
# E266: arithmetic operator without whitespace on both sides


.PHONY: install test format lint docs checktodo