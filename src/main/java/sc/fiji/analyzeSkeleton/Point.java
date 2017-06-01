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

/**
 * This class represents a 3D point or position on a 3D image. Therefore, the 
 * coordinates are integer.
 */
public class Point 
{
	/** x- coordinate */
	public int x = 0;
	/** y- coordinate */
	public int y = 0;
	/** z- coordinate */
	public int z = 0;

	/**
	 * Create point from integer coordinates.
	 * 
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate
	 */
	public Point(int x, int y, int z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/**
	 * Convert point to string.
	 */
	public String toString(){
		return "(" + this.x + ", " + this.y + ", " + this.z + ")";
	}

	/**
	 * Override equals method to compare points.
	 * @param o input object
	 * @return true if the input object is equal to this Point
	 */
	public boolean equals(Object o)
	{
		if (this == o) return true;

		if (o == null || getClass() != o.getClass()) return false;


		final Point p = (Point) o;
		return p.x == this.x && p.y == this.y && p.z == this.z;
	}

	@Override
    public Point clone() {
	    return new Point(x, y, z);
    }
}// end class point
