
SHELL := /bin/bash

# Points to Utility Directory
COMMON_REPO = ./
ABS_COMMON_REPO = $(shell readlink -f $(COMMON_REPO))

UTILS_PATH = ./utils


include $(UTILS_PATH)/help.mk
include $(UTILS_PATH)/utils.mk

#export  XCL_EMULATION_MODE=sw_emu
TARGETS := hw
TARGET  := $(TARGETS)
DEVICES := xilinx_vcu1525_xdma_201830_1
# device list:
# xilinx_vcu1525_xdma_201830_1
# xilinx_u200_xdma_201830_2
# xilinx_u250_xdma_201830_2

DEVICE  := $(DEVICES)


app := pr
APP = $(app)

APPCONFIG = ./application/$(APP)


include $(APPCONFIG)/config.mk

include $(APPCONFIG)/build.mk

include ./application/common.mk
CXX=gcc


.PHONY:application
application:: $(APPCONFIG)
	@touch $(APPCONFIG)/build.mk
	@echo $(opencl_LDFLAGS)

.PHONY: all clean cleanall docs emconfig
all: $(EXECUTABLE) $(BINARY_CONTAINERS) emconfig application

.PHONY: exe
exe: cleanexe approximation_motifs_scheme_1  approximation_triangle_scheme_1  approximation_motifs_scheme_3 approximation_triangle_scheme_3

# Building kernel

include ./application/common_gs_kernel.mk

ifeq ($(strip $(HAVE_APPLY)), $(strip $(VAR_TRUE)))
$(XCLBIN)/vertexApply.$(TARGET).$(DSA).xo: $(APPLY_KERNEL_PATH)/vertex_apply.cpp
	mkdir -p $(XCLBIN)
	$(XOCC) $(CLFLAGS) -c -k vertexApply -I'$(<D)' -o'$@' '$<'
include ./application/common_apply_kernel.mk
include $(APPCONFIG)/apply_kernel.mk
endif


$(XCLBIN)/graph_fpga.$(TARGET).$(DSA).xclbin: $(BINARY_CONTAINER_OBJS)
	$(XOCC) $(CLFLAGS) -l $(LDCLFLAGS)   $(BINARY_LINK_OBJS) -o'$@' $(+)

# Building Host
$(EXECUTABLE): $(HOST_SRCS)
	mkdir -p $(XCLBIN)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -o '$@' $(LDFLAGS)

.PHONY: approximation_motifs_scheme_3
approximation_motifs_scheme_3: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_motifs_scheme_3 -DMAX_RUN_STEPS=12 -DSUB_EST=1000  -o '$@' $(LDFLAGS)


.PHONY: approximation_triangle_scheme_3
approximation_triangle_scheme_3: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_triangle_scheme_3 -DMAX_RUN_STEPS=12 -DSUB_EST=1000 -o '$@' $(LDFLAGS)


.PHONY: approximation_triangle_scheme_4
approximation_triangle_scheme_4: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_triangle_scheme_4 -DMAX_RUN_STEPS=21 -DSUB_EST=1 -o '$@' $(LDFLAGS)



.PHONY: approximation_motifs_scheme_1
approximation_motifs_scheme_1: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_motifs_scheme_1 -DSUB_EST=1 -o '$@' $(LDFLAGS)

.PHONY: approximation_motifs_scheme_2
approximation_motifs_scheme_2: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_motifs_scheme_2 -DSUB_EST=1 -o '$@' $(LDFLAGS)

.PHONY: approximation_triangle_scheme_1
approximation_triangle_scheme_1: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_triangle_scheme_1 -DSUB_EST=1 -o '$@' $(LDFLAGS)

.PHONY: approximation_triangle_scheme_2
approximation_triangle_scheme_2: $(HOST_SRCS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -DAPPROXIMATE_FUNCTION=approximation_triangle_scheme_2 -DSUB_EST=1 -o '$@' $(LDFLAGS)


emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(DEVICE) --od $(EMCONFIG_DIR)

.PHONY: hwemuprepare
hwemuprepare:
ifeq ($(TARGET),$(filter $(TARGET), hw_emu))
	@echo "prepare for hw_emu"
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	$(CP) $(UTILS_PATH)/sdaccel.ini .
	source $(UTILS_PATH)/hw_emu.sh
else
	@echo "prepare for hw"
endif


check: all 

ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) ./$(EXECUTABLE) 
else
	 ./$(EXECUTABLE)
endif
	sdx_analyze profile -i sdaccel_profile_summary.csv -f html

# Cleaning stuff



cleanexe:
	-$(RMDIR) approximation_motifs_*
	-$(RMDIR) approximation_triangle_*
clean:
	-$(RMDIR) $(EXECUTABLE) $(XCLBIN)/{*sw_emu*,*hw_emu*} 
	-$(RMDIR) sdaccel_* TempConfig system_estimate.xtxt *.rpt
	-$(RMDIR) src/*.ll _xocc_* .Xil emconfig.json dltmp* xmltmp* *.log *.jou *.wcfg *.wdb
	-$(RMDIR) .Xil
	-$(RMDIR) *.zip
	-$(RMDIR) *.str

cleanall: clean
	-$(RMDIR) $(XCLBIN)
	-$(RMDIR) ./_x
	-$(RMDIR) ./membership.out


cleandir: cleanall
	-$(RMDIR) host_graph_fpga*
	-$(RMDIR) xclbin*
	-$(RMDIR) .run