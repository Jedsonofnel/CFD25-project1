.PHONY: all bin/cli

all: bin/cli

bin/cli:
	go build -o ./bin/cli .
