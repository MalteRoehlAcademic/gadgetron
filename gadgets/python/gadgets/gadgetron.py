try:
    import GadgetronPythonMRI
    import ismrmrd
except ImportError:
    pass

import time
import numpy as np

# Gateway for python to call back c++ or next gadget in python
class Gadget(object):
    def __init__(self, next_gadget=None):
        self.next_gadget = next_gadget
        self.params = dict()
        self.results = []

    def set_parameter(self, name, value):
        self.params[name] = value

    def get_parameter(self, name):
        return self.params.get(name, None)

    def __call__(self, *args):
        self.process(args)
        return self.get_results()

    def set_next_gadget(self, gadget):
        self.next_gadget = gadget

    def process_config(self, conf):
        pass

    def process(self, header, *args):
        # do work here
        self.put_next(header,*args)

    def wait(self):
        pass
    
    def put_next(self, *args):
        if self.next_gadget is not None:
            if isinstance(self.next_gadget, Gadget):
                """
                if len(args) == 3 and not isinstance(args[2],ismrmrd.Meta): #Data with meta data we assume
                    meta = ismrmrd.Meta()
                    meta = ismrmrd.Meta.deserialize(args[2])
                    new_args = (args[0], args[1], meta)
                    self.next_gadget.process(*new_args)
                elif len(args) == 3 and not isinstance(args[2],ismrmrd.Meta): #Data with meta data we assume
                    meta = ismrmrd.Meta()
                    meta = ismrmrd.Meta.deserialize(args[2])
                    new_args = (args[0], args[1], meta)
                    self.next_gadget.process(*new_args)
                else:
                """
                self.next_gadget.process(*args)
            elif isinstance(self.next_gadget, GadgetronPythonMRI.GadgetReference):
                if len(args) > 3:
                    raise Exception("Only two or 3 return arguments are currently supported when returning to Gadgetron framework")
                if isinstance(args[0], ismrmrd.AcquisitionHeader):
                    if len(args) == 3 and (isinstance(args[2], ismrmrd.Meta) == False):
                        self.next_gadget.return_acquisition_traj(args[0],args[1].astype('complex64'),args[2].astype('float32'))
                    else:
                        self.next_gadget.return_acquisition(args[0],args[1].astype('complex64'))
                elif isinstance(args[0],IsmrmrdImageArray): 
                    self.next_gadget.return_ismrmrd_image_array(args[0])
                elif isinstance(args[0], ismrmrd.ImageHeader):
                    header = args[0]
                    if (args[1].dtype == np.uint16):
                        if len(args) == 3:
                            if isinstance(args[2], ismrmrd.Meta):
                                self.next_gadget.return_image_ushort_attr(header,args[1], args[2].serialize())
                            else:
                                self.next_gadget.return_image_ushort_attr(header,args[1], str(args[2]))
                        else:
                            self.next_gadget.return_image_ushort(header,args[1])
                    elif (args[1].dtype == np.float32):
                        if len(args) == 3:
                            if isinstance(args[2], ismrmrd.Meta):
                                self.next_gadget.return_image_float_attr(header, args[1], args[2].serialize())
                            else:
                                self.next_gadget.return_image_float_attr(header, args[1], str(args[2]))
                        else:
                            self.next_gadget.return_image_float(header,args[1])
                    else:
                        if len(args) == 3:
                            if isinstance(args[2], ismrmrd.Meta):
                                self.next_gadget.return_image_cplx_attr(header, args[1].astype('complex64'), str(args[2].serialize()))
                            else:
                                self.next_gadget.return_image_cplx_attr(header, args[1].astype('complex64'), str(args[2]))
                        else:
                            self.next_gadget.return_image_cplx(header,args[1].astype('complex64'))
                elif len(args[0]) > 0 and isinstance(args[0][0],IsmrmrdReconBit): 
                    self.next_gadget.return_recondata(args[0])
                elif len(args[0]) > 0 and isinstance(args[0][0],GenericSpiralReconJob): 
                    self.next_gadget.return_recondataspiral(args[0])
                else:
                    raise("Unsupported types when returning to Gadgetron framework")
            else:
                raise("next_gadget is set to unsupported type")
        else:
            self.results.append(list(args))

    def get_results(self):
        results = self.results
        self.results = []
        return results


class WrappedGadget(object):
    def __init__(self, dllname, classname, gadgetname):
        self.gadgetname = gadgetname
        self.classname = classname
        self.dllname = dllname
        self.parameters = dict()
    
class WrapperGadget(Gadget):
    
    def __init__(self, dllname, classname, gadgetname=None, next_gadget=None):
        if gadgetname == None:
            gadgetname = classname
        Gadget.__init__(self, next_gadget)
        self.controller_ = GadgetronPythonMRI.GadgetInstrumentationStreamController()
        self.controller_.prepend_gadget(gadgetname,dllname,classname)
        self.controller_.set_python_gadget(self)
        self.wrapped_gadgets = list()
        self.wrapped_gadgets.append(WrappedGadget(dllname,classname,gadgetname))

    def prepend_gadget(self,dllname, classname, gadgetname=None):
        self.controller_.prepend_gadget(gadgetname,dllname,classname)
        self.wrapped_gadgets.insert(0,WrappedGadget(dllname,classname,gadgetname))

    def wait(self):
       self.controller_.close()

    def process_config(self, conf):
        self.controller_.put_config(conf)
        return 0

    def process(self, header, *args):
        if len(args) > 2:
            raise("Only two or three arguments are currently supported when sending data to Gadgetron framework")
        if isinstance(header, ismrmrd.AcquisitionHeader):
            if len(args) == 2 and isinstance(args[1], ismrmrd.Meta):
                self.controller_.put_acquisition(header,args[0].astype('complex64'))
            else:
                self.controller_.put_acquisition_traj(header,args[0].astype('complex64'),args[1].astype('float32'))
        elif isinstance(header, ismrmrd.ImageHeader):
            if (args[0].dtype == np.uint16):
                if len(args) == 2:
                    self.controller_.put_image_ushort_attr(header,args[0], args[1].serialize())
                else:
                    self.controller_.put_image_ushort(header,args[0])
            elif (args[0].dtype == np.float32):
                if len(args) == 2:
                    self.controller_.put_image_float_attr(header, args[0], args[1].serialize())
                else:
                    self.controller_.put_image_float(header,args[0])
            else:   
                if len(args) == 2:
                    self.controller_.put_image_cplx_attr(header, args[0].astype('complex64'), args[1].serialize())
                else:
                    self.controller_.put_image_cplx(header,args[0].astype('complex64'))
        elif hasattr(args[0],"__getitem__") and hasattr(args[0][0],"headers"):
            self.controller_.put_ismrmrd_image_array(args[0])
        elif hasattr(args[0],"__getitem__") and hasattr(args[0][0],"data"):
            
            if hasattr(args[0],"__getitem__") and hasattr(args[0][0],"density"):
                self.controller_.put_recondataspiral(args[0])
            else:
                self.controller_.put_recondata(args[0])
        else:
            raise("Unsupported types when sending data to Gadgetron framework")
        return 0

    def set_parameter(self,gadgetname,parameter,value):
        for g in self.wrapped_gadgets:
            if g.gadgetname == gadgetname:
                g.parameters[parameter] = value
                break

        self.controller_.set_parameter(gadgetname,parameter,value)

class FunctionGadget(Gadget):
    """A Gadget with a configurable `process` function.

    Params:
        fn: `process` function
    """
    def __init__(self, fn, next_gadget=None):
        super(FunctionGadget, self).__init__(next_gadget)
        self.process = fn


def gadget_chain_wait(first_gadget):
    g = first_gadget;
    while (g):
        g.wait()
        g = g.next_gadget

def gadget_chain_config(first_gadget, conf):
    g = first_gadget;
    while (g):
        g.process_config(conf)
        g = g.next_gadget
        
def get_last_gadget(first_gadget):
    g = first_gadget;
    while (True):
        if g.next_gadget:
            g = g.next_gadget
        else:
            break
    return g

class SamplingLimit:
    def __init__(self):
        self.min = 0
        self.center = 0
        self.max = 0

class SpiralSamplingLimit:
    def __init__(self):
        self.min = 0
        self.center = 0
        self.max = 0

class SegmentedMapping:
    def __init__(self):
        self.min = 0
        self.max = 0
        self.lower_boundary_bin = 0
        self.upper_boundary_bin = 0

class SamplingDescription:
    def __init__(self):
        self.encoded_FOV = (0.0,0.0,0.0)
        self.recon_FOV = (0.0,0.0,0.0)
        self.encoded_matrix = (0,0,0)
        self.recon_matrix = (0,0,0)
        self.sampling_limits = (SamplingLimit(),SamplingLimit(),SamplingLimit())

class SpiralSamplingDescription:
    def __init__(self):
        self.encoded_FOV = (0,0,0)
        self.recon_FOV = (0,0,0)
        self.encoded_resolution = 0.0
        self.recon_resolution = 0.0
        self.flags = 0
        self.interleaves = 0
        self.oversampling_factor = 0
        self.kernel_width = 0
        self.matrix_size_factor = 0
        self.adc_sampling_time = 0
        self.acceleration = 0
        self.channels = 0
        self.samples_per_interleave = 0
        self.matrix = (0,0,0)
        self.segmented_mapping = (SegmentedMapping(),SegmentedMapping())
        self.total_number_of_segments = (0,0)
        self.current_number_segments = (0,0)
        self.sampling_limits = SpiralSamplingLimit()

        self.slices = ()
        self.averages = ()
        self.repetitions = ()
 


class IsmrmrdDataBuffered:
    def __init__(self,data,headers,sampling=SamplingDescription(),trajectory=None):
        self.data = data 
        if (trajectory is not None): 
            self.trajectory =trajectory 
        self.headers =headers 
        self.sampling = sampling

class GenericSpiralReconJob:
    def __init__(self,data,headers,trajectory,density,sampling=SpiralSamplingDescription(),result=None,field_map=None,t2_star_map=None,csm=None,mask=None,motion_field=None,inverse_motion_field=None,reg=None):
        self.data = data
        self.headers = headers
        self.trajectory = trajectory
        self.density = density
        self.sampling = sampling
        self.result = result
        self.field_map = field_map
        self.t2_star_map = t2_star_map
        self.csm = csm
        self.mask = mask
        self.motion_field = motion_field
        self.inverse_motion_field = inverse_motion_field
        self.reg = reg

        self.number_of_slices = len(sampling.slices)
        self.number_of_averages = len(sampling.averages)
        self.number_of_repetitions = len(sampling.repetitions)


        if result is None:
            self.has_result_ = False
        else:
            self.has_result_ = True

        if field_map is None:
            self.has_field_map_ = False
        else:
            self.has_field_map_ = True

        if t2_star_map is None:
            self.has_t2_star_map_ = False
        else:
            self.has_t2_star_map_ = True

        if csm is None:
            self.has_csm_ = False
        else:
            self.has_csm_ = True

        if reg is None:
            self.has_reg_ = False
        else:
            self.has_reg_ = True

        if motion_field is None:
            self.has_motion_field_ = False
        else:
            self.has_motion_field_ = True

        if inverse_motion_field is None:
            self.has_inverse_motion_field_ = False
        else:
            self.has_inverse_motion_field_ = True

        if mask is None:
            self.has_mask_ = False
        else:
            self.has_mask_ = True

    


        
        

class IsmrmrdReconBit:
    def __init__(self,data,ref=None):
        self.data = data
        if (ref != None): 
            self.ref = ref

class IsmrmrdImageArray: 
    def __init__(self,data=None,headers=None,meta=None,waveform=None,acq_headers=None):
        self.data=data
        self.headers=headers
        self.meta=meta
        self.waveform=waveform
        self.acq_headers=acq_headers
