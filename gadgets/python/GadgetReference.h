#pragma once

#include "Gadget.h"
#include "gadgetronpython_export.h"

#include <ismrmrd/ismrmrd.h>
#include <boost/python.hpp>
#include "mri_core_data.h"
#include "Generic_Spiral_Recon_Job.h"

namespace Gadgetron{

  class EXPORTGADGETSPYTHON GadgetReference
  {

  public:
    GadgetReference();
    ~GadgetReference();

    int set_gadget(Gadget* g)
    {
      gadget_ = g;
      return 0;
    }

    template<class TH, class TD> int return_data(TH header, boost::python::object arr, const char* meta = 0);
    template<class TH, class TD, class TE> int return_data_acq(TH header, boost::python::object arr, boost::python::object arr2);
    int return_acquisition(ISMRMRD::AcquisitionHeader acq, boost::python::object arr);
    int return_acquisition_traj(ISMRMRD::AcquisitionHeader acq, boost::python::object arr, boost::python::object arr2);
    int return_recondata(boost::python::object arr);
    int return_recondataspiral(boost::python::object arr);
    int return_ismrmrd_image_array(boost::python::object rec);
    int return_image_cplx(ISMRMRD::ImageHeader img, boost::python::object arr);
    int return_image_cplx_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);
    int return_image_float(ISMRMRD::ImageHeader img, boost::python::object arr);
    int return_image_float_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);
    int return_image_ushort(ISMRMRD::ImageHeader img, boost::python::object arr);
    int return_image_ushort_attr(ISMRMRD::ImageHeader img, boost::python::object arr, const char* meta);

  protected:
    Gadget* gadget_;
  };
}
