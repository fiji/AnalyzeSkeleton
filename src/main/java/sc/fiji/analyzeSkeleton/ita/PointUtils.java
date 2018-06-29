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

import java.util.Collection;

import org.joml.Vector3d;

import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Utility methods for the classes of the ita package.
 *
 * @author Richard Domander (Royal Veterinary College, London)
 * @author Alessandro Felder (Royal Veterinary College, London)
 */
public final class PointUtils {

	private PointUtils() {}

	/**
	 * Returns the center of the given points.
	 *
	 * @param points points of a {@link Vertex}.
	 * @return {x, y, z} coordinates of the centroid.
	 */
	public static Vector3d centroid(final Collection<Point> points) {
		final Vector3d centroid = new Vector3d();
		points.forEach(p -> centroid.add(p.x, p.y, p.z));
		centroid.div(points.size());
		return centroid;
	}
}
