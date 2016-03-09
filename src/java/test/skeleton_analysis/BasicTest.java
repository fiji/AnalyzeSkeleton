package skeleton_analysis;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

public class BasicTest {

	/**
	 * Main method to test and debug the AnalyzeSkeleton plugin
	 *  
	 * @param args
	 */
	public static void main( final String[] args )
	{
		ImageJ.main( args );
		
		ImagePlus imp = IJ.openImage( 
				AnalyzeSkeleton_.class.getResource( 
						"/bat-cochlea-skeleton.zip" ).getFile() );
		imp.show();
		AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();		
		skel.setup("", imp);		
		skel.run( null );
	}
	
}
