/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Wed Oct 23 18:08:19 2019
 * Arguments: 
 * "-B""macro_default""-W""java:generate_sending_signal,Test""-T""link:lib""-d""D:\\gitRepo\\OFDM\\generate_sending_signal\\for_testing""class{Test:D:\\gitRepo\\OFDM\\decode_received_signal.m,D:\\gitRepo\\OFDM\\generate_sending_signal.m}"
 */

package generate_sending_signal;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class Generate_sending_signalMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "generate_sen_1E598847500A983B9F87778C9C9BB2AC";
    
    /** Component name */
    private static final String sComponentName = "generate_sending_signal";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(Generate_sending_signalMCRFactory.class)
        );
    
    
    private Generate_sending_signalMCRFactory()
    {
        // Never called.
    }
    
    public static MWMCR newInstance(MWComponentOptions componentOptions) throws MWException
    {
        if (null == componentOptions.getCtfSource()) {
            componentOptions = new MWComponentOptions(componentOptions);
            componentOptions.setCtfSource(sDefaultComponentOptions.getCtfSource());
        }
        return MWMCR.newInstance(
            componentOptions, 
            Generate_sending_signalMCRFactory.class, 
            sComponentName, 
            sComponentId,
            new int[]{9,5,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
