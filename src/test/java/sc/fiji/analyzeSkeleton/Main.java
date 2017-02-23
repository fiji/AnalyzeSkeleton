/*
 * #%L
 * AnalyzeSkeleton_ plugin for ImageJ.
 * %%
 * Copyright (C) 2008 - 2017 Ignacio Arganda-Carreras.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package sc.fiji.analyzeSkeleton;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

public class Main {

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
