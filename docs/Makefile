API_DIR='api'

all: configure parse plot

configure: api_dir
	pydocmd simple hydrad_tools.configure++ hydrad_tools.configure.Configure++ > $(API_DIR)/configure.md

parse: api_dir
	pydocmd simple hydrad_tools.parse++ hydrad_tools.parse.Strand++ hydrad_tools.parse.Profile++ > $(API_DIR)/parse.md

plot: api_dir
	pydocmd simple hydrad_tools.visualize+ hydrad_tools.visualize.plot_strand hydrad_tools.visualize.animate_strand > $(API_DIR)/visualize.md

api_dir:
	mkdir -p $(API_DIR)

clean:
	echo "Cleaning up API docs"
	rm -r $(API_DIR)/*