BSC_DIR = basics/
THIRD_PARTY_DIR = $(BSC_DIR)third_party/
IMGUI_DIR = $(THIRD_PARTY_DIR)imgui/

fetregister: dear_imgui
		make -C fetregister

fetbenchmark:
		make -C fetbenchmark

dear_imgui:
		make -C $(IMGUI_DIR)
