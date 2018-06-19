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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.joml.Vector3d;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Utility methods for analysing the junctions of a {@link Graph}.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
public final class VertexUtils {

	private VertexUtils() {}

	/**
	 * Calculates the angles between the branches of each vertex.
	 * <p>
	 * For example, if a vertex has 3 branches, it's returned as a list of 3
	 * angles. That is, angles between branches 1-2, 1-3 and 2-3. A vertex with no
	 * branches or only one branch returns an empty list.
	 * </p>
	 * <p>
	 * Angle between branches 1 and 2 is measured from their respective end points
	 * (vertices opposite the current "junction", i.e. the vertex being
	 * processed).
	 * </p>
	 *
	 * @param vertices vertices from a graph.
	 * @return a list of angles in radians of the vertex branches.
	 */
	public static List<List<Double>> getNJunctionAngles(
		final Collection<Vertex> vertices)
	{
		final List<List<Double>> graphAngles = new ArrayList<>(vertices.size());
		for (final Vertex vertex : vertices) {
			final List<Double> junctionAngles = getAngles(vertex);
			graphAngles.add(junctionAngles);
		}
		return graphAngles;
	}

	/**
	 * Groups vertices by their valence, i.e. by their number of branches.
	 * <p>
	 * Only vertices with valence in the given range are grouped.
	 * </p>
	 *
	 * @param vertices vertices from a graph
	 * @param min minimum valance (inclusive) for grouped vertices
	 * @param max maximum valence (inclusive) for grouped vertices
	 * @return a (Valence, Vertices with valence) mapping.
	 * @throws IllegalArgumentException if min &lt; 0, or min &gt max.
	 */
	public static Map<Integer, List<Vertex>> groupByValence(
		final Collection<Vertex> vertices, final int min, final int max)
		throws IllegalArgumentException
	{
		if (min < 0) {
			throw new IllegalArgumentException("Minimum must be non-negative");
		}
		if (min > max) {
			throw new IllegalArgumentException(
				"Minimum must be less or equal to maximum");
		}
		final Predicate<Vertex> isInValenceRange = v -> {
			final int valence = v.getBranches().size();
			return min <= valence && valence <= max;
		};
		return vertices.stream().filter(isInValenceRange).collect(Collectors
			.groupingBy(v -> v.getBranches().size()));
	}

	// region -- Helper methods --
	private static List<Double> getAngles(final Vertex junction) {
		final List<Double> angles = new ArrayList<>();
		final Vector3d centroid = Util.centroid(junction.getPoints());
		final List<Edge> branches = junction.getBranches();
		for (int i = 0; i < branches.size() - 1; i++) {
			for (int j = i + 1; j < branches.size(); j++) {
				final Vector3d endpoint = getOppositeCentroid(junction, centroid,
					branches.get(i));
				final Vector3d endpoint2 = getOppositeCentroid(junction, centroid,
					branches.get(j));
				final double angle = endpoint.angle(endpoint2);
				angles.add(angle);
			}
		}
		return angles;
	}

	private static Vector3d getOppositeCentroid(final Vertex vertex,
		final Vector3d centroid, final Edge edge)
	{
		final List<Point> points = edge.getOppositeVertex(vertex).getPoints();
		final Vector3d endPoint = Util.centroid(points);
		return endPoint.sub(centroid);
	}
	// endregion
}
