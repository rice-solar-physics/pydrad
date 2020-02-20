API_DIR='api'

all: configure parse plot

configure: api_dir
	pydocmd simple pydrad.configure++ pydrad.configure.Configure++ > $(API_DIR)/configure.md

parse: api_dir
	pydocmd simple pydrad.parse++ pydrad.parse.Strand++ pydrad.parse.Profile++ > $(API_DIR)/parse.md

plot: api_dir
	pydocmd simple pydrad.visualize+ pydrad.visualize.plot_strand pydrad.visualize.animate_strand pydrad.visualize.plot_profile pydrad.visualize.plot_time_distance > $(API_DIR)/visualize.md

api_dir:
	mkdir -p $(API_DIR)

clean:
	echo "Cleaning up API docs"
	rm -r $(API_DIR)/*