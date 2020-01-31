#pragma once
#include "python_toolbox.h"
#include "python_numpy_wrappers.h"

#include "hoNDArray.h"
#include "Generic_Spiral_Recon_Job.h"
#include "log.h"

#include <boost/python.hpp>
namespace bp = boost::python;

namespace Gadgetron {

class IsmrmrdReconDataSpiral_to_python_object {
public:
  static PyObject* convert(const IsmrmrdReconDataSpiral & reconDataSpiral) {
//      initialize_python();
      GILLock lock;
    bp::object pygadgetron = bp::import("gadgetron");


    
    auto pyReconData = bp::list();
    for (auto & reconBit : reconDataSpiral.rbit_ )
    {
      auto data = SpiralReconDataToPython(reconBit);
      pyReconData.append(data);
    }

    // increment the reference count so it exists after `return`
    return bp::incref(pyReconData.ptr());
  }

private:
  static bp::object SpiralReconDataToPython( const GenericSpiralReconJob & dataSpiral)
  {
    bp::object pygadgetron = bp::import("gadgetron");

    auto data = bp::object(dataSpiral.data_);
    auto headers = boost::python::object(dataSpiral.headers_);
    auto trajectory = bp::object(dataSpiral.trajectory_);
    auto density = bp::object(dataSpiral.density_);
    auto result = bp::object(dataSpiral.result_) ;
    auto field_map =  bp::object(dataSpiral.field_map_) ;
    auto t2_star_map =  bp::object(dataSpiral.t2_star_map_) ;
    auto csm =  bp::object(dataSpiral.csm_) ;
    auto mask = bp::object(dataSpiral.mask_) ;
    auto motion_field = bp::object(dataSpiral.motion_field_) ;
    auto inverse_motion_field = bp::object(dataSpiral.inverse_motion_field_) ;
    auto reg = bp::object(dataSpiral.reg_) ;
    auto sampling = SpiralSamplingDescriptionToPython(dataSpiral.sampling_);

    auto buffer = pygadgetron.attr("GenericSpiralReconJob")(data,headers,trajectory,density,sampling,result,field_map,t2_star_map,csm,mask
    ,motion_field,inverse_motion_field,reg);

    bp::incref(data.ptr());
    bp::incref(headers.ptr());
    bp::incref(trajectory.ptr());
    bp::incref(density.ptr());
    bp::incref(result.ptr());
    bp::incref(field_map.ptr());
    bp::incref(t2_star_map.ptr());
    bp::incref(csm.ptr());
    bp::incref(mask.ptr());
    bp::incref(motion_field.ptr());
    bp::incref(inverse_motion_field.ptr());
    bp::incref(reg.ptr());
    bp::incref(sampling.ptr());

    return buffer;
  }

  static bp::object SpiralSamplingDescriptionToPython(const SpiralSamplingDescription & sD)
  {
    bp::object pygadgetron = bp::import("gadgetron");
    bp::object result;
    try 
    {
      auto tmp = pygadgetron.attr("SpiralSamplingDescription");
      result = tmp();
      result.attr("encoded_FOV") = bp::make_tuple(sD.encoded_FOV_[0],sD.encoded_FOV_[1],sD.encoded_FOV_[2]);
      result.attr("recon_FOV") = bp::make_tuple(sD.recon_FOV_[0],sD.recon_FOV_[1],sD.recon_FOV_[2]);
      result.attr("matrix") = bp::make_tuple(sD.matrix_[0],sD.matrix_[1],sD.matrix_[2]);

      result.attr("encoded_resolution") = sD.encoded_resolution_;
      result.attr("recon_resolution") = sD.recon_resolution_;
      result.attr("flags") = sD.flags_;
      result.attr("interleaves") = sD.interleaves_;
      result.attr("oversampling_factor") = sD.oversampling_factor_;
      result.attr("kernel_width") = sD.kernel_width_;
      result.attr("matrix_size_factor") = sD.matrix_size_factor_;
      result.attr("adc_sampling_time") = sD.adc_sampling_time_;
      result.attr("acceleration") = sD.acceleration_;
      result.attr("channels") = sD.channels_;
      result.attr("samples_per_interleave") = sD.samples_per_interleave_;
      result.attr("total_number_of_segments") = bp::make_tuple(sD.total_number_of_segments_[0],sD.total_number_of_segments_[1]);
      result.attr("current_number_segments") = bp::make_tuple(sD.current_number_segments_[0],sD.current_number_segments_[1]);

    } 
    catch (bp::error_already_set const &)
    {
        std::string err = pyerr_to_string();
        GERROR(err.c_str());
        throw std::runtime_error(err);
    }


    auto sl = bp::list();
    for (size_t ii =0; ii < sD.slices_.size();ii++)
      sl.append(sD.slices_[ii]);

    auto av = bp::list();
    for (size_t ii =0; ii < sD.averages_.size();ii++)
      av.append(sD.averages_[ii]);

    auto re = bp::list();
    for (size_t ii =0; ii < sD.repetitions_.size();ii++)
      re.append(sD.repetitions_[ii]);


    result.attr("slices") = sl;
    result.attr("averages") = av;
    result.attr("repetitions") = re;

    std::vector<bp::object> SLe;
    for (int i = 0; i < 2; i++)
    {
      auto segmented_mapping = pygadgetron.attr("SegmentedMapping")();

      segmented_mapping.attr("min") = sD.segmented_mapping_[i].min_;
      segmented_mapping.attr("max") = sD.segmented_mapping_[i].max_;
      segmented_mapping.attr("lower_boundary_bin") = sD.segmented_mapping_[i].lower_boundary_bin_;
      segmented_mapping.attr("upper_boundary_bin") = sD.segmented_mapping_[i].upper_boundary_bin_;
      SLe.push_back(segmented_mapping);
    }
    result.attr("segmented_mapping") = bp::make_tuple(SLe[0],SLe[1]);

    auto sampling_limit = pygadgetron.attr("SpiralSamplingLimit")();
    sampling_limit.attr("min") = sD.sampling_limits_.min_;
    sampling_limit.attr("max") = sD.sampling_limits_.max_;
    sampling_limit.attr("center") = sD.sampling_limits_.center_;
  
    result.attr("sampling_limit") = sampling_limit;

    return result;
  }
};


/// Used for making an hoNDArray from a NumPy array
struct IsmrmrdReconDataSpiral_from_python_object {
  IsmrmrdReconDataSpiral_from_python_object() {
    // actually register this converter with Boost
    bp::converter::registry::push_back(
        &convertible,
        &construct,
        bp::type_id<IsmrmrdReconDataSpiral>());
  }

  /// Returns NULL if the object is not convertible. Or well.... it should
  static void* convertible(PyObject* obj) {
    return obj;
  }

  /// Construct an hoNDArray in-place
  static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
    void* storage = ((bp::converter::rvalue_from_python_storage<IsmrmrdReconData >*)data)->storage.bytes;
    IsmrmrdReconDataSpiral* reconData = new (storage) IsmrmrdReconDataSpiral;
    data->convertible = storage;


    try {
      bp::list pyRecondata((bp::handle<>(bp::borrowed(obj))));
      auto length = bp::len(pyRecondata);
      GDEBUG("Recon data spiral length: %i\n",length);
      for (int i = 0; i < length; i++){
        bp::object reconBit = pyRecondata[i];
        GenericSpiralReconJob rBit = extractSpiralJob(reconBit);
        reconData->rbit_.push_back(rBit);
      }

    }catch (const bp::error_already_set&) {
      std::string err = pyerr_to_string();
      GERROR(err.c_str());
      throw std::runtime_error(err);
    }
  }

  static GenericSpiralReconJob extractSpiralJob(bp::object pyDataSpiral)
  {
    GenericSpiralReconJob result;

    result.headers_ = bp::extract<hoNDArray<ISMRMRD::AcquisitionHeader>>(pyDataSpiral.attr("headers"));
    result.data_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataSpiral.attr("data"));
    result.trajectory_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("trajectory"));
    result.density_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("density"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"result"))
      result.result_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataSpiral.attr("result"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"field_map"))
      result.field_map_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("field_map"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"t2_star_map"))
      result.t2_star_map_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("t2_star_map"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"csm"))
      result.csm_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataSpiral.attr("csm"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"reg"))
      result.reg_ = bp::extract<hoNDArray<std::complex<float>>>(pyDataSpiral.attr("reg"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"motion_field"))
      result.motion_field_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("motion_field"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"inverse_motion_field"))
      result.inverse_motion_field_ = bp::extract<hoNDArray<float>>(pyDataSpiral.attr("inverse_motion_field"));

    if (PyObject_HasAttrString(pyDataSpiral.ptr(),"mask"))
      result.mask_ = bp::extract<hoNDArray<unsigned short>>(pyDataSpiral.attr("mask"));

    auto pySampling = pyDataSpiral.attr("sampling");
    SpiralSamplingDescription sampling;


    for (int i = 0; i < 3; i++)
      sampling.encoded_FOV_[i] = bp::extract<float>(pySampling.attr("encoded_FOV")[i]);
    for (int i = 0; i < 3; i++)
      sampling.matrix_[i] = bp::extract<uint16_t>(pySampling.attr("matrix")[i]);
    for (int i = 0; i < 3; i++)
      sampling.recon_FOV_[i] = bp::extract<float>(pySampling.attr("recon_FOV")[i]);

    sampling.encoded_resolution_ = bp::extract<float>(pySampling.attr("encoded_resolution"));
    sampling.recon_resolution_ = bp::extract<float>(pySampling.attr("recon_resolution"));
    sampling.flags_ = bp::extract<int32_t>(pySampling.attr("flags"));
    sampling.interleaves_ = bp::extract<uint16_t>(pySampling.attr("interleaves"));
    sampling.oversampling_factor_ = bp::extract<float>(pySampling.attr("oversampling_factor"));
    sampling.kernel_width_ = bp::extract<float>(pySampling.attr("kernel_width"));
    sampling.matrix_size_factor_ = bp::extract<float>(pySampling.attr("matrix_size_factor"));
    sampling.adc_sampling_time_ = bp::extract<float>(pySampling.attr("adc_sampling_time"));
    sampling.acceleration_ = bp::extract<uint16_t>(pySampling.attr("acceleration"));
    sampling.channels_ = bp::extract<uint16_t>(pySampling.attr("channels"));
    sampling.samples_per_interleave_ = bp::extract<uint16_t>(pySampling.attr("samples_per_interleave"));

    for (int i = 0; i < 2; i++)
      sampling.total_number_of_segments_[i] = bp::extract<uint16_t>(pySampling.attr("total_number_of_segments")[i]);
    for (int i = 0; i < 2; i++)
      sampling.current_number_segments_[i] = bp::extract<uint16_t>(pySampling.attr("current_number_segments")[i]);



    auto pySlices = pyDataSpiral.attr("slices");
    auto pyAverages = pyDataSpiral.attr("averages");
    auto pyRepetitions = pyDataSpiral.attr("repetitions");

    auto lengthsl = bp::len(pySlices);
    auto lengthav = bp::len(pyAverages);
    auto lengthrep = bp::len(pyRepetitions);

    sampling.averages_.resize(lengthav);
    sampling.slices_.resize(lengthsl);
    sampling.repetitions_.resize(lengthrep);

    auto sl = bp::list();
    for (size_t ii =0; ii < lengthsl;ii++)
      sampling.slices_[ii] = bp::extract<uint16_t>(pySampling.attr("slices")[ii]);

    auto av = bp::list();
    for (size_t ii =0; ii < lengthav;ii++)
      sampling.averages_[ii] = bp::extract<uint16_t>(pySampling.attr("averages")[ii]);

    auto re = bp::list();
    for (size_t ii =0; ii < lengthrep;ii++)
      sampling.repetitions_[ii] = bp::extract<uint16_t>(pySampling.attr("repetitions")[ii]);



    auto pySLs = pySampling.attr("sampling_limit");
    sampling.sampling_limits_.min_ = bp::extract<uint16_t>(pySLs.attr("min"));
    sampling.sampling_limits_.center_ = bp::extract<uint16_t>(pySLs.attr("center"));
    sampling.sampling_limits_.max_ = bp::extract<uint16_t>(pySLs.attr("max"));

    for (int i = 0; i < 2; i++){
      auto pySL = pySampling.attr("segmented_mapping")[i];
      sampling.segmented_mapping_[i].min_ = bp::extract<float>(pySL.attr("min"));
      sampling.segmented_mapping_[i].max_ = bp::extract<float>(pySL.attr("max"));
      sampling.segmented_mapping_[i].lower_boundary_bin_ = bp::extract<float>(pySL.attr("lower_boundary_bin"));
      sampling.segmented_mapping_[i].upper_boundary_bin_ = bp::extract<float>(pySL.attr("upper_boundary_bin"));
    }

    result.sampling_ = sampling;
    return result;
  }


};


/// Partial specialization of `python_converter` for hoNDArray
template<> struct python_converter<IsmrmrdReconDataSpiral> {
  static void create()
  {
    // register hoNDArray converter
    bp::type_info info = bp::type_id<IsmrmrdReconDataSpiral >();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    // only register if not already registered!
    if (nullptr == reg || nullptr == (*reg).m_to_python) {
      bp::to_python_converter<IsmrmrdReconDataSpiral, IsmrmrdReconDataSpiral_to_python_object >();
      IsmrmrdReconDataSpiral_from_python_object();
    }

  }
};

}
