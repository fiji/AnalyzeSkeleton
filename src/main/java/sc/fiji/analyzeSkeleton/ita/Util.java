/*
 * #%L
 * Graph analysis tools for the AnalyzeSkeleton_ plugin for ImageJ.
 * %%
 * Copyright (C) 2017 - 2018 Richard Domander, Alessandro Felder & Michael Doube
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
package sc.fiji.analyzeSkeleton.ita;

import org.joml.Vector3d;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;

import java.util.Collection;

/**
 * Utility methods for the classes of the ita package.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
final class Util {
	private Util() {}

	/**
	 * Returns the center of the given points.
	 *
	 * @param points points of vertices in a {@link Graph}.
	 * @return {x, y, z} coordinates of the centroid.
	 */
	static Vector3d centroid(final Collection<Point> points) {
		final Vector3d centroid = new Vector3d();
		points.forEach(p -> centroid.add(p.x, p.y, p.z));
		centroid.div(points.size());
		return centroid;
	}
}
