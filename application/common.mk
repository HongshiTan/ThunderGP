
VAR_TRUE=true

CP = cp -rf

XCLBIN := ./xclbin_$(APP)
DSA := $(call device2sandsa, $(DEVICE))

CXX := $(XILINX_SDX)/bin/xcpp
XOCC := $(XILINX_SDX)/bin/xocc



GS_KERNEL_PATH    = ./libfpga/common
ifeq ($(strip $(CUSTOMIZE_APPLY)), $(strip $(VAR_TRUE)))
APPLY_KERNEL_PATH = $(APPCONFIG)
else
APPLY_KERNEL_PATH = ./libfpga/common
endif


include $(ABS_COMMON_REPO)/utils/opencl.mk



#--kernel_frequency <frequency>
#--xp prop:solution.kernel_compiler_margin=<Frequency Percentage>
#--xp param:compiler.enableAutoFrequencyScaling=0

HOST_SRCS = ./host_graph.cpp ./libgraph/graph.cpp ./libgraph/data_helper.cpp

ifeq ($(strip $(GP_SYSTEM)), $(strip $(VAR_TRUE)))

HOST_SRCS += ./libgraph/memory/he_mem.cpp 

HOST_SRCS += ./libgraph/host_graph_sw_mem.cpp
HOST_SRCS += ./libgraph/host_graph_sw_partition.cpp
HOST_SRCS += ./libgraph/host_graph_sw_kernel.cpp

HOST_SRCS += ./libgraph/host_graph_sw_dataflow.cpp

HOST_SRCS += ./libgraph/scheduler/host_graph_scheduler.cpp
HOST_SRCS += ./libgraph/scheduler/$(SCHEDULER)/scheduler.cpp

HOST_SRCS += $(APPCONFIG)/l3DataPrepare.cpp


ifeq ($(strip $(HAVE_APPLY)), $(strip $(VAR_TRUE)))

ifeq ($(strip $(CUSTOMIZE_APPLY)), $(strip $(VAR_TRUE)))
HOST_SRCS += $(APPCONFIG)/host_vertex_apply.cpp
endif

HOST_SRCS += ./libgraph/verification/host_graph_verification_apply.cpp
endif

HOST_SRCS += ./libgraph/verification/host_graph_verification_gs.cpp

endif

# Host compiler global settings
CXXFLAGS := $(opencl_CXXFLAGS) -Wall
CXXFLAGS += -I/$(XILINX_SDX)/Vivado_HLS/include/ -O3 -g -fmessage-length=0 -std=c++14 -Wno-deprecated-declarations
CXXFLAGS += -I ./
CXXFLAGS += -I ./libfpga
CXXFLAGS += -I ./libgraph
CXXFLAGS += -I ./libgraph/memory
CXXFLAGS += -I ./libgraph/scheduler
CXXFLAGS += -I ./libgraph/verification
CXXFLAGS += -I $(APPCONFIG)

# Host linker flags
LDFLAGS := $(opencl_LDFLAGS)
LDFLAGS += -lrt -lstdc++  -lxilinxopencl
LDFLAGS += -lgmp -lboost_system -lboost_random

ifeq ($(TARGET),$(filter $(TARGET), hw_emu))
CLFLAGS = -g -t $(TARGET)
else
CLFLAGS = -t $(TARGET)
endif

# Kernel compiler global settings

CLFLAGS += --platform $(DEVICE) --save-temps  -O3
CLFLAGS += -I ./
CLFLAGS += -I ./libfpga
CLFLAGS += -I $(APPCONFIG)
CLFLAGS += --xp prop:solution.kernel_compiler_margin=10%

# Kernel linker flags
LDCLFLAGS += --xp prop:solution.kernel_compiler_margin=10% --kernel_frequency=$(FREQ)

EXECUTABLE = host_graph_fpga_$(APP)

EMCONFIG_DIR = $(XCLBIN)/$(DSA)

BINARY_CONTAINERS += $(XCLBIN)/graph_fpga.$(TARGET).$(DSA).xclbin


#Include Libraries

include $(UTILS_PATH)/xcl/xcl.mk
CXXFLAGS +=  $(xcl_CXXFLAGS)
LDFLAGS +=   $(xcl_CXXFLAGS)
HOST_SRCS += $(xcl_SRCS)


#############################################################################
#                                                                           #
#                     Specific Build Configuration                          #
#                                                                           #
#############################################################################

ifeq ($(strip $(HAVE_FULL_SLR)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DSUB_PARTITION_NUM=4
else
CXXFLAGS += -DSUB_PARTITION_NUM=1
endif

ifeq ($(strip $(HAVE_APPLY)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DHAVE_APPLY=1
else
CXXFLAGS += -DHAVE_APPLY=0
endif


ifeq ($(strip $(HAVE_VERTEX_ACTIVE_BIT)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DHAVE_VERTEX_ACTIVE_BIT=1
CLFLAGS  += -DHAVE_VERTEX_ACTIVE_BIT=1
else
CXXFLAGS += -DHAVE_VERTEX_ACTIVE_BIT=0
CLFLAGS  += -DHAVE_VERTEX_ACTIVE_BIT=0
endif

ifeq ($(strip $(HAVE_EDGE_PROP)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DHAVE_EDGE_PROP=1
CLFLAGS  += -DHAVE_EDGE_PROP=1
else
CXXFLAGS += -DHAVE_EDGE_PROP=0
CLFLAGS  += -DHAVE_EDGE_PROP=0
endif


ifeq ($(strip $(HAVE_UNSIGNED_PROP)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DHAVE_UNSIGNED_PROP=1
CLFLAGS  += -DHAVE_UNSIGNED_PROP=1
else
CXXFLAGS += -DHAVE_UNSIGNED_PROP=0
CLFLAGS  += -DHAVE_UNSIGNED_PROP=0
endif


ifeq ($(strip $(HAVE_APPLY_OUTDEG)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DHAVE_APPLY_OUTDEG=1
CLFLAGS  += -DHAVE_APPLY_OUTDEG=1
else
CXXFLAGS += -DHAVE_APPLY_OUTDEG=0
CLFLAGS  += -DHAVE_APPLY_OUTDEG=0
endif



ifeq ($(strip $(CUSTOMIZE_APPLY)), $(strip $(VAR_TRUE)))
CXXFLAGS += -DCUSTOMIZE_APPLY=1
CLFLAGS  += -DCUSTOMIZE_APPLY=1
else
CXXFLAGS += -DCUSTOMIZE_APPLY=0
CLFLAGS  += -DCUSTOMIZE_APPLY=0
endif


#############################################################################
#                                                                           #
#                     Specific Configuration                                #
#                                                                           #
#############################################################################

CLFLAGS  += -DQUEUE_SIZE_FILTER=$(QUEUE_SIZE_FILTER)
CLFLAGS  += -DQUEUE_SIZE_MEMORY=$(QUEUE_SIZE_MEMORY)
CLFLAGS  += -DLOG_SCATTER_CACHE_BURST_SIZE=$(LOG_SCATTER_CACHE_BURST_SIZE)