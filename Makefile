.PHONY: all bin/cli peclet-sweep grid-sweep

all: bin/cli

bin/cli:
	go build -o ./bin/cli .

peclet-sweep: bin/cli
	@echo "Running Peclet number sweep (Pe = 0.1 to 50)..."
	@mkdir -p data
	@for u in 0.1 0.2 0.5 1 2 5 10 20 50; do \
		echo "  Pe = $$u"; \
		./bin/cli -n 20 -u $$u > data/pe$$u.csv; \
	done
	@echo "Done! Data saved to data/"

grid-sweep: bin/cli
	@echo "Running grid refinement sweep (n = 5 to 160)..."
	@mkdir -p data
	@for n in 5 10 20 40 80 160; do \
		echo "  n = $$n"; \
		./bin/cli -n $$n -u 2.5 > data/grid$$n.csv; \
	done
	@echo "Done! Data saved to data/"

clean:
	rm -rf bin/ data/
