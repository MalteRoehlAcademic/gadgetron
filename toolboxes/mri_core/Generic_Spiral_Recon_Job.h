#ifndef GenericSpiralReconJob_H
#define GenericSpiralReconJob_H

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/waveform.h"
#include "ismrmrd/meta.h"
#include <vector>
#include <set>
#include "hoNDArray.h"
#include <boost/optional.hpp>
#include <complex>


namespace Gadgetron
{
  
  class SpiralSamplingLimit
  {
  public:
    int16_t min_;
    int16_t center_;
    int16_t max_;

    SpiralSamplingLimit()
    {
        min_ = 0;
        center_ = 0;
        max_ = 0;
    }
  };

  class SegmentedMapping
  {
  public:
    float min_;
    float max_;

    float lower_boundary_bin_;
    float upper_boundary_bin_;

    SegmentedMapping()
    {
        min_ = 0;
        max_ = 0;
        lower_boundary_bin_ = 0;
        upper_boundary_bin_ = 0;
    }
  };
  
  class SpiralSamplingDescription
  {
  public:
    // encoding FOV
    float encoded_FOV_[3];
    // recon FOV
    float recon_FOV_[3];

    float encoded_resolution_;
    float recon_resolution_;

    int32_t flags_;

    uint16_t interleaves_;
    std::vector<uint16_t> slices_;
    std::vector<uint16_t> averages_;
    std::vector<uint16_t> repetitions_;

    float oversampling_factor_;
    float kernel_width_;
    float matrix_size_factor_;

    float adc_sampling_time_;

    uint16_t acceleration_;
    uint16_t channels_;
    uint16_t samples_per_interleave_;

    uint16_t matrix_[3];

    SegmentedMapping segmented_mapping_[2];
    uint16_t total_number_of_segments_[2];
    uint16_t current_number_segments_[2];
    
    // sampled range along Readout, Interleaves
    // min, max and center
    SpiralSamplingLimit sampling_limits_;

    SpiralSamplingDescription()
    {
        encoded_FOV_[0] = 0;
        encoded_FOV_[1] = 0;
        encoded_FOV_[2] = 0;

        recon_FOV_[0] = 0;
        recon_FOV_[1] = 0;
        recon_FOV_[2] = 0;

        flags_ = 0;
        adc_sampling_time_ = 0;
        encoded_resolution_ = 0;
        recon_resolution_ = 0;
        oversampling_factor_ = 0;
        kernel_width_ = 0;
        acceleration_ = 0;
        matrix_size_factor_ = 0;
        samples_per_interleave_ = 0;
        interleaves_ = 0;
        slices_.clear();
        averages_.clear();
        repetitions_.clear();
        channels_ = 0;
        total_number_of_segments_[0] = 0;
        total_number_of_segments_[1] = 0;
        current_number_segments_[0] = 0;
        current_number_segments_[1] = 0;

        matrix_[0] = 0;
        matrix_[1] = 0;
        matrix_[2] = 0;
    }

    SpiralSamplingDescription(SpiralSamplingDescription& obj)
    {
        encoded_FOV_[0] = obj.encoded_FOV_[0];
        encoded_FOV_[1] = obj.encoded_FOV_[1];
        encoded_FOV_[2] = obj.encoded_FOV_[2];

        recon_FOV_[0] = obj.recon_FOV_[0];
        recon_FOV_[1] = obj.recon_FOV_[1];
        recon_FOV_[2] = obj.recon_FOV_[2];

        flags_ = obj.flags_;
        adc_sampling_time_ = obj.adc_sampling_time_;
        encoded_resolution_ = obj.encoded_resolution_;
        recon_resolution_ = obj.recon_resolution_;
        oversampling_factor_ = obj.oversampling_factor_;
        kernel_width_ = obj.kernel_width_;
        acceleration_ = obj.acceleration_;
        segmented_mapping_[0] = obj.segmented_mapping_[0];
        segmented_mapping_[1] = obj.segmented_mapping_[1];
        matrix_size_factor_ = obj.matrix_size_factor_;
        samples_per_interleave_ = obj.samples_per_interleave_;
        interleaves_ = obj.interleaves_;
        slices_ = obj.slices_;
        averages_ = obj.averages_;
        repetitions_ = obj.repetitions_;
        channels_ = obj.channels_;
        total_number_of_segments_[0] = obj.total_number_of_segments_[0] ;
        total_number_of_segments_[1] = obj.total_number_of_segments_[1] ;
        current_number_segments_[0] = obj.current_number_segments_[0];
        current_number_segments_[1] = obj.current_number_segments_[1];

        matrix_[0] = obj.matrix_[0];
        matrix_[1] = obj.matrix_[1];
        matrix_[2] = obj.matrix_[2];
    }

    ~SpiralSamplingDescription() {}
  };

  struct GenericSpiralReconJob
  {
  public:

    //4D, fixed order [Interleaves, Average, SLICE, Repetition]
    hoNDArray< ISMRMRD::AcquisitionHeader > headers_;

    //7D, fixed order [Samples, Channels, Interleaves, Average, SLICE, Repetition]
    hoNDArray< std::complex<float> > data_;
    
    //7D, fixed order [X,Y,Samples, Interleaves, Average, SLICE, Repetition]
    hoNDArray<float> trajectory_;

    // 6D, density weights [Samples, Interleaves, Average, SLICE, Repetition]
    hoNDArray<float>  density_;

    //7D, fixed order [X, Y, Average, SLICE, Repetition]
    hoNDArray<std::complex<float>> result_;

    //7D, fixed order [X, Y, SLICE]
    hoNDArray<float> field_map_;

    //7D, fixed order [X, Y, SLICE]
    hoNDArray<float> t2_star_map_;

    //7D, fixed order [X, Y,CHA, SLICE]
    hoNDArray<std::complex<float>> csm_;

    //7D, fixed order [X, Y, Average, SLICE, Repetition]
    hoNDArray<std::complex<float>> reg_;

    //7D, fixed order [X, Y, Average, SLICE, Repetition]
    hoNDArray<float> motion_field_;

    //7D, fixed order [X, Y, Average, SLICE, Repetition]
    hoNDArray<float> inverse_motion_field_;

    //7D, fixed order [X, Y, SLICE]
    hoNDArray<float> mask_;

    bool has_result_ = false;
    bool has_field_map_ = false;
    bool has_t2_star_map_ = false;
    bool has_csm_ = false;
    bool has_reg_ = false;
    bool has_motion_field_ = false;
    bool has_inverse_motion_field_ = false;
    bool has_mask_ = false;

    SpiralSamplingDescription sampling_;

    GenericSpiralReconJob() {}
    GenericSpiralReconJob(const GenericSpiralReconJob& obj)
    {

        if (obj.data_.get_data_ptr())
        {
            if (this->data_.get_data_ptr())
                if (this->data_.delete_data_on_destruct()) this->data_.clear();
             else 
                this->data_ = hoNDArray<std::complex<float>>();
            this->data_.copyFrom(obj.data_);
        }

        if (obj.trajectory_.get_data_ptr())
        {
            if (this->trajectory_.get_data_ptr())
                if (this->trajectory_.delete_data_on_destruct()) this->trajectory_.clear();
             else 
                this->trajectory_ = hoNDArray<float>();
            this->trajectory_.copyFrom(obj.trajectory_);
        }

        if (obj.density_.get_data_ptr())
        {
            if (this->density_.get_data_ptr())
                if (this->density_.delete_data_on_destruct()) this->density_.clear();
             else 
                this->density_ = hoNDArray<float>();
            this->density_.copyFrom(obj.density_);
        }

        if (obj.field_map_.get_data_ptr())
        {
            if (this->field_map_.get_data_ptr())
                if (this->field_map_.delete_data_on_destruct()) this->field_map_.clear();
             else 
                this->field_map_ = hoNDArray<float>();
            this->field_map_.copyFrom(obj.field_map_);
        }

        if (obj.t2_star_map_.get_data_ptr())
        {
            if (this->t2_star_map_.get_data_ptr())
                if (this->t2_star_map_.delete_data_on_destruct()) this->t2_star_map_.clear();
             else 
                this->t2_star_map_ = hoNDArray<float>();
            this->t2_star_map_.copyFrom(obj.t2_star_map_);
        }

        if (obj.csm_.get_data_ptr())
        {
            if (this->csm_.get_data_ptr())
                if (this->csm_.delete_data_on_destruct()) this->csm_.clear();
             else 
                this->csm_ = hoNDArray<std::complex<float>>();
            this->csm_.copyFrom(obj.csm_);
        }


                
        if (obj.reg_.get_data_ptr())
        {
            if (this->reg_.get_data_ptr())
                if (this->reg_.delete_data_on_destruct()) this->reg_.clear();
             else 
                this->reg_ = hoNDArray<std::complex<float>>();
            this->reg_.copyFrom(obj.reg_);
        }

        if (obj.motion_field_.get_data_ptr())
        {
            if (this->motion_field_.get_data_ptr())
                if (this->motion_field_.delete_data_on_destruct()) this->motion_field_.clear();
             else 
                this->motion_field_ = hoNDArray<float>();
            this->motion_field_.copyFrom(obj.motion_field_);
        }

        if (obj.inverse_motion_field_.get_data_ptr())
        {
            if (this->inverse_motion_field_.get_data_ptr())
                if (this->inverse_motion_field_.delete_data_on_destruct()) this->inverse_motion_field_.clear();
             else 
                this->inverse_motion_field_ = hoNDArray<float>();
            this->inverse_motion_field_.copyFrom(obj.inverse_motion_field_);
        }
        if (obj.result_.get_data_ptr())
        {
            if (this->result_.get_data_ptr())
                if (this->result_.delete_data_on_destruct()) this->result_.clear();
             else 
                this->result_ = hoNDArray<std::complex<float>>();
            this->result_.copyFrom(obj.result_);
        }

        if (obj.mask_.get_data_ptr())
        {
            if (this->mask_.get_data_ptr())
                if (this->mask_.delete_data_on_destruct()) this->mask_.clear();
             else 
                this->mask_ = hoNDArray<float>();
            this->mask_.copyFrom(obj.mask_);
        }   

        
        this->headers_.copyFrom(obj.headers_);
        this->sampling_ = obj.sampling_;
    }

    GenericSpiralReconJob& operator=(GenericSpiralReconJob other)
    {

        if (other.data_.get_data_ptr())
        {
            if (this->data_.get_data_ptr())
                if (this->data_.delete_data_on_destruct()) this->data_.clear();
             else 
                this->data_ = hoNDArray<std::complex<float>>();
            this->data_.copyFrom(other.data_);
        }
        else
        {
            this->data_ = hoNDArray<std::complex<float>>();
        }
        
        if (other.trajectory_.get_data_ptr())
        {
            if (this->trajectory_.get_data_ptr())
                if (this->trajectory_.delete_data_on_destruct()) this->trajectory_.clear();
             else 
                this->trajectory_ = hoNDArray<float>();
            this->trajectory_.copyFrom(other.trajectory_);
        }
        else
        {
            this->trajectory_ = hoNDArray<float>();
        }

        if (other.density_.get_data_ptr())
        {
            if (this->density_.get_data_ptr())
                if (this->density_.delete_data_on_destruct()) this->density_.clear();
             else 
                this->density_ = hoNDArray<float>();
            this->density_.copyFrom(other.density_);
        }
        else
        {
            this->density_ = hoNDArray<float>();
        }

        if (other.field_map_.get_data_ptr())
        {
            if (this->field_map_.get_data_ptr())
                if (this->field_map_.delete_data_on_destruct()) this->field_map_.clear();
             else 
                this->field_map_ = hoNDArray<float>();
            this->field_map_.copyFrom(other.field_map_);
        }
        else
        {
            this->field_map_ = hoNDArray<float>();
        }

        if (other.t2_star_map_.get_data_ptr())
        {
            if (this->t2_star_map_.get_data_ptr())
                if (this->t2_star_map_.delete_data_on_destruct()) this->t2_star_map_.clear();
             else 
                this->t2_star_map_ = hoNDArray<float>();
            this->t2_star_map_.copyFrom(other.t2_star_map_);
        }
        else
        {
            this->t2_star_map_ = hoNDArray<float>();
        }

        if (other.csm_.get_data_ptr())
        {
            if (this->csm_.get_data_ptr())
                if (this->csm_.delete_data_on_destruct()) this->csm_.clear();
             else 
                this->csm_ = hoNDArray<std::complex<float>>();
            this->csm_.copyFrom(other.csm_);
        }
        else
        {
            this->csm_ = hoNDArray<std::complex<float>>();
        }

                
        if (other.reg_.get_data_ptr())
        {
            if (this->reg_.get_data_ptr())
                if (this->reg_.delete_data_on_destruct()) this->reg_.clear();
             else 
                this->reg_ = hoNDArray<std::complex<float>>();
            this->reg_.copyFrom(other.reg_);
        }
        else
        {
            this->reg_ = hoNDArray<std::complex<float>>();
        }

        if (other.motion_field_.get_data_ptr())
        {
            if (this->motion_field_.get_data_ptr())
                if (this->motion_field_.delete_data_on_destruct()) this->motion_field_.clear();
             else 
                this->motion_field_ = hoNDArray<float>();
            this->motion_field_.copyFrom(other.motion_field_);
        }
        else
        {
            this->motion_field_ = hoNDArray<float>();
        }

        if (other.inverse_motion_field_.get_data_ptr())
        {
            if (this->inverse_motion_field_.get_data_ptr())
                if (this->inverse_motion_field_.delete_data_on_destruct()) this->inverse_motion_field_.clear();
             else 
                this->inverse_motion_field_ = hoNDArray<float>();
            this->inverse_motion_field_.copyFrom(other.inverse_motion_field_);
        }        
        else
        {
            this->inverse_motion_field_ = hoNDArray<float>();
        }

        if (other.result_.get_data_ptr())
        {
            if (this->result_.get_data_ptr())
                if (this->result_.delete_data_on_destruct()) this->result_.clear();
             else 
                this->result_ = hoNDArray<std::complex<float>>();
            this->result_.copyFrom(other.result_);
        }
        else
        {
            this->result_ = hoNDArray<std::complex<float>>();
        }

        if (other.mask_.get_data_ptr())
        {
            if (this->mask_.get_data_ptr())
                if (this->mask_.delete_data_on_destruct()) this->mask_.clear();
             else 
                this->mask_ = hoNDArray<float>();
            this->mask_.copyFrom(other.mask_);
        }
        else
        {
            this->mask_ = hoNDArray<float>();
        }

        
        this->headers_.copyFrom(other.headers_);
        this->sampling_ = other.sampling_;

        return *this;
    }

    ~GenericSpiralReconJob() {}

    void clear()
    {

        if (this->data_.get_data_ptr())
        {
            if (this->data_.delete_data_on_destruct()) this->data_.clear();
        }
        if (this->trajectory_.get_data_ptr())
        {
            if (this->trajectory_.delete_data_on_destruct()) this->trajectory_.clear();
        }

        if (this->density_.get_data_ptr())
        {
            if (this->density_.delete_data_on_destruct()) this->density_.clear();
        }

        if (this->field_map_.get_data_ptr())
        {
            if (this->field_map_.delete_data_on_destruct()) this->field_map_.clear();
        }

        if (this->t2_star_map_.get_data_ptr())
        {
            if (this->t2_star_map_.delete_data_on_destruct()) this->t2_star_map_.clear();
        }

        if (this->csm_.get_data_ptr())
        {
            if (this->csm_.delete_data_on_destruct()) this->csm_.clear();
        }

   

        if (this->reg_.get_data_ptr())
        {
            if (this->reg_.delete_data_on_destruct()) this->reg_.clear();
        }

        if (this->motion_field_.get_data_ptr())
        {
            if (this->motion_field_.delete_data_on_destruct()) this->motion_field_.clear();
        }

        if (this->inverse_motion_field_.get_data_ptr())
        {
            if (this->inverse_motion_field_.delete_data_on_destruct()) this->inverse_motion_field_.clear();
        }

        if (this->result_.get_data_ptr())
        {
            if (this->result_.delete_data_on_destruct()) this->result_.clear();
        }

        if (this->mask_.get_data_ptr())
        {
            if (this->mask_.delete_data_on_destruct()) this->mask_.clear();
        }


        if (this->headers_.delete_data_on_destruct()) headers_.clear();
    }
  };

  struct IsmrmrdReconDataSpiral
  {
    public:
        std::vector<GenericSpiralReconJob> rbit_;
        IsmrmrdReconDataSpiral() {}
        IsmrmrdReconDataSpiral(const IsmrmrdReconDataSpiral& obj)
        {
            if (!this->rbit_.empty())
                this->rbit_.clear();
            this->rbit_ = obj.rbit_;
        }
        ~IsmrmrdReconDataSpiral() {}
  };

}
#endif

