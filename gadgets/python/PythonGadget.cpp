#include "PythonGadget.h"

namespace Gadgetron {

    int PythonGadget::process(ACE_Message_Block* mb)
    {
        GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* hma = AsContainerMessage<ISMRMRD::AcquisitionHeader>(mb);
        if (hma)
        {
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* dmb = AsContainerMessage< hoNDArray< std::complex<float> > >(hma->cont());
            if (!dmb)
            {
                GERROR("Error converting data array from message block for Acquisition\n");
                hma->release();
                return GADGET_FAIL;;
            }
            GadgetContainerMessage< hoNDArray<float> >* traj = AsContainerMessage< hoNDArray< float > >(hma->cont()->cont());
            if (traj)
            {
                return this->process(hma, dmb, traj);
            }
            else
            {
                GadgetContainerMessage< ISMRMRD::MetaContainer>* mmb = AsContainerMessage< ISMRMRD::MetaContainer >(dmb->cont());
                if (!mmb)
                    mmb = new  GadgetContainerMessage< ISMRMRD::MetaContainer>();
                return this->process(hma, dmb, mmb);
            }

        }
        else
        {
            GadgetContainerMessage<ISMRMRD::ImageHeader>* hmi = AsContainerMessage<ISMRMRD::ImageHeader>(mb);
            if (hmi)
                return this->process_image(hmi);
        }

        {
            auto recon_data = AsContainerMessage<IsmrmrdReconData>(mb);
            if (recon_data) {
                GDEBUG("Calling into python gadget with IsmrmrdReconData\n");
                return this->process(recon_data);
            }
        }

        {
            auto recon_data = AsContainerMessage<IsmrmrdReconDataSpiral>(mb);
            if (recon_data) {
                GDEBUG("Calling into python gadget with IsmrmrdReconDataSpiral\n");
                return this->process(recon_data);
            }
        }

        {
            auto array_data = AsContainerMessage<IsmrmrdImageArray>(mb);
            if (array_data)
            {
                GDEBUG("Calling into python gadget with IsmrmrdImageArray\n");
                return this->process(array_data);
            }
        }

        // since most python gadgets do not handle waveform yet 
        // and by default, waveform is passed down the chain (this is special), 
        // these lines can cause failure of python recon
        // TODO: a better solution for python waveform handling
        {
            if (auto waveform_header = AsContainerMessage<ISMRMRD::ISMRMRD_WaveformHeader>(mb)){
                if (auto waveform_data = AsContainerMessage<hoNDArray<uint32_t>>(mb->cont())){
                    return this->process(waveform_header,waveform_data);
                }

            }
        }

        if (pass_on_undesired_data.value()) {
            return this->next()->putq(mb);
        }
        else {
            GERROR("%s. This is neither an acquisition or an image. Something is wrong here", class_name.c_str());

            mb->release();
            return GADGET_FAIL;
        }
    }

    int PythonGadget::process_image(GadgetContainerMessage<ISMRMRD::ImageHeader>* hmi)
    {
        ISMRMRD::ImageHeader* h = hmi->getObjectPtr();
        GadgetContainerMessage< ISMRMRD::MetaContainer>* mmb = 0;

        if (hmi->cont()) {
            mmb = AsContainerMessage< ISMRMRD::MetaContainer >(hmi->cont()->cont());
            if (!mmb)
                mmb = new  GadgetContainerMessage< ISMRMRD::MetaContainer>();
        }

        switch (h->data_type) {
        case (ISMRMRD::ISMRMRD_USHORT):
            return this->process(hmi, AsContainerMessage< hoNDArray< uint16_t > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_SHORT):
            return this->process(hmi, AsContainerMessage< hoNDArray< int16_t > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_UINT):
            return this->process(hmi, AsContainerMessage< hoNDArray< uint32_t > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_INT):
            return this->process(hmi, AsContainerMessage< hoNDArray< int32_t > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_FLOAT):
            return this->process(hmi, AsContainerMessage< hoNDArray< float > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_DOUBLE):
            return this->process(hmi, AsContainerMessage< hoNDArray< double > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_CXFLOAT):
            return this->process(hmi, AsContainerMessage< hoNDArray< std::complex<float> > >(hmi->cont()), mmb);
            break;
        case (ISMRMRD::ISMRMRD_CXDOUBLE):
            return this->process(hmi, AsContainerMessage< hoNDArray< std::complex<double> > >(hmi->cont()), mmb);
            break;
        default:
            GERROR("Unknown image data_type %d received\n", h->data_type);
            hmi->release();
            return GADGET_FAIL;
            break;
        }
    }
    GADGET_FACTORY_DECLARE(PythonGadget)
}
