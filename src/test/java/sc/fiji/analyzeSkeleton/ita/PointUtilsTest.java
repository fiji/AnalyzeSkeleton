
package sc.fiji.analyzeSkeleton.ita;

import static java.util.Arrays.asList;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.joml.Vector3d;
import org.junit.Test;

import sc.fiji.analyzeSkeleton.Point;

/**
 * Tests for {@link PointUtils}.
 *
 * @author Richard Domander (Royal Veterinary College, London)
 * @author Alessandro Felder (Royal Veterinary College, London)
 */
public class PointUtilsTest {

	@Test
	public void testCentroid() {
		final Vector3d expected = new Vector3d(0.5, 0.5, 0.5);
		final List<Point> unitCube = asList(new Point(0, 0, 0), new Point(1, 0, 0),
			new Point(1, 1, 0), new Point(0, 1, 0), new Point(0, 0, 1), new Point(1,
				0, 1), new Point(1, 1, 1), new Point(0, 1, 1));

		final Vector3d result = PointUtils.centroid(unitCube);

		assertEquals(expected, result);
	}

	@Test
	public void testCentroidEmptyList() {
		final Vector3d expected = new Vector3d(Double.NaN, Double.NaN, Double.NaN);
		final List<Point> empty = new ArrayList<>();

		final Vector3d result = PointUtils.centroid(empty);

		assertEquals(expected, result);
	}

	@Test
	public void testCentroidOnePoint() {
		final Vector3d expected = new Vector3d(1.0, 2.0, 3.0);
		final List<Point> points = Collections.singletonList(new Point(1, 2, 3));

		final Vector3d result = PointUtils.centroid(points);

		assertEquals(expected, result);
	}
}
