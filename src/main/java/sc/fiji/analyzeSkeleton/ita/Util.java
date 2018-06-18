package sc.fiji.analyzeSkeleton.ita;

import org.joml.Vector3d;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;

import java.util.Collection;

/**
 * @author Richard Domander
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
