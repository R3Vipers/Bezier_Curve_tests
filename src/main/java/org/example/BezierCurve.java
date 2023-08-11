package org.example;

import java.awt.geom.Point2D;
import java.util.ArrayList;


public class BezierCurve {

    // Calculate the binomial coefficient (n choose k)
    private static int binomialCoefficient(int n, int k) {
        if (k == 0 || k == n) {
            return 1;
        } else {
            return binomialCoefficient(n - 1, k - 1) + binomialCoefficient(n - 1, k);
        }
    }

    // Calculate a single Bernstein polynomial
    private static double bernsteinPolynomial(int i, int n, double t) {
        return binomialCoefficient(n, i) * Math.pow(1 - t, n - i) * Math.pow(t, i);
    }

    // Generate points on a Bézier curve
    public static Point2D.Double[] generateBezierCurve(Point2D.Double[] controlPoints, int numPoints) {
        int n = controlPoints.length - 1;
        Point2D.Double[] curvePoints = new Point2D.Double[numPoints];

        for (int j = 0; j < numPoints; j++) {
            double t = (double) j / (numPoints - 1);
            double x = 0.0;
            double y = 0.0;

            for (int i = 0; i <= n; i++) {
                double b = bernsteinPolynomial(i, n, t);
                x += controlPoints[i].getX() * b;
                y += controlPoints[i].getY() * b;
            }

            curvePoints[j] = new Point2D.Double(x, y);
        }

        return curvePoints;
    }

    public static ArrayList <Vector2d> generateVelocityVector(Point2D.Double[] curvePoints){
        ArrayList <Vector2d> velocityVectors = new ArrayList<Vector2d>();
        for(int i = 1; i < curvePoints.length-2; i++) {
            double dx = (curvePoints[i-1].getX()-curvePoints[i].getX());
            double dy = (curvePoints[i-1].getY()-curvePoints[i].getY());
            double dt = 1.0/curvePoints.length;
            double dxdt = dx/dt;
            double dydt = dy/dt;
            Vector2d vector = new Vector2d(dxdt, dydt);
            velocityVectors.add(vector);
        }

        return velocityVectors;
    }

    public static double [] generateCentripetalForceVectorMagnitude(Point2D.Double[] curvePoints, double mass, double velocity) {
        double [] centripetalForceVectorMagnitude = new double[curvePoints.length];
        for(int i = 0; i < curvePoints.length-2; i++) {
            double dx1 = (curvePoints[i+1].getX()-curvePoints[i].getX());
            double dx2 = (curvePoints[i+2].getX()-curvePoints[i+1].getX());
            double dy1 = (curvePoints[i+1].getY()-curvePoints[i].getY());
            double dy2 = (curvePoints[i+2].getY()-curvePoints[i+1].getY());
            double size = curvePoints.length;
            double dt = 1.0/size;
            double dxdt1 = dx1/dt;
            double dydt1 = dy1/dt;
            double dxdt2 = dx2/dt;
            double dydt2 = dy2/dt;
            double d2xdt2 = (dxdt2-dxdt1)/dt;
            double d2ydt2 = (dydt2-dydt1)/dt;
            double dsec = (((dxdt1*d2ydt2)-(dydt1*d2xdt2))/Math.pow(dxdt1, 2.0))/dxdt1;
            double mper = -dxdt1/dydt1;
            double r = -(Math.pow(Math.pow(1 + Math.pow(mper, 2.0), 0.5)/mper,3.0)/dsec);
            //double r = (Math.pow((Math.pow(dxdt2, 2.0) + Math.pow(dydt2, 2.0)), 3/2)) / ((dxdt2*d2ydt2)-(dydt2*d2xdt2));
            double force = mass*((Math.pow(velocity, 2))/r);
            centripetalForceVectorMagnitude [i] = r;
        }
        return centripetalForceVectorMagnitude;
    }

    public static Point2D.Double[] generateCircleCenter(Point2D.Double[] curvePoints) {
        Point2D.Double[] circleCenter = new Point2D.Double[curvePoints.length-1];
        double [] radius = new double[curvePoints.length-1];
        for(int i = 0; i < curvePoints.length-2; i++) {
            double dx1 = (curvePoints[i+1].getX()-curvePoints[i].getX());
            double dx2 = (curvePoints[i+2].getX()-curvePoints[i+1].getX());
            double dy1 = (curvePoints[i+1].getY()-curvePoints[i].getY());
            double dy2 = (curvePoints[i+2].getY()-curvePoints[i+1].getY());
            double size = curvePoints.length;
            double dt = 1.0/size;
            double dxdt1 = dx1/dt;
            double dydt1 = dy1/dt;
            double dxdt2 = dx2/dt;
            double dydt2 = dy2/dt;
            double d2xdt2 = (dxdt2-dxdt1)/dt;
            double d2ydt2 = (dydt2-dydt1)/dt;
            double dsec = (((dxdt1*d2ydt2)-(dydt1*d2xdt2))/Math.pow(dxdt1, 2.0))/dxdt1;
            double mper = -dxdt1/dydt1;
            double r = -(Math.pow(Math.pow(1 + Math.pow(mper, 2.0), 0.5)/mper,3.0)/dsec);
            radius [i] = r;
        }

        for(int i = 0; i < curvePoints.length-1; i++) {
            double dx = (curvePoints[i+1].getX()-curvePoints[i].getX());
            double dy = (curvePoints[i+1].getY()-curvePoints[i].getY());
            double dt = 1.0/curvePoints.length;
            double dxdt = dx/dt;
            double dydt = dy/dt;
            double dxtdyt = -dxdt/dydt;
            double x = curvePoints[i].getX()- (radius[i] * Math.cos(Math.atan(dxtdyt)));
            double y = curvePoints[i].getY()+ (radius[i] * Math.sin(-Math.atan(dxtdyt)));
            circleCenter[i] = new Point2D.Double(x, y);
        }
        return  circleCenter;
    }

    // Example usage
    public static void main(String[] args) {
        Point2D.Double[] controlPoints = {
                new Point2D.Double(0.0, 0.0),
                new Point2D.Double(5.00, 10.00),
                new Point2D.Double(10.00, -5.00),
                new Point2D.Double(15.00, 0.0)
        };

        int numPoints = 1000;
        double mass = 11.3398;
        double velocity = 1.2;
        Point2D.Double[] curvePoints = generateBezierCurve(controlPoints, numPoints);
        ArrayList velocityVectors = generateVelocityVector(curvePoints);
        double [] velocityMagnitude = generateCentripetalForceVectorMagnitude(curvePoints, mass, velocity);
        Point2D.Double[] circleCenters = generateCircleCenter(curvePoints);

        // Print the generated points
        System.out.println("points: \n");
        for (int i = 0; i < numPoints; i++) {
            System.out.println("Point " + i + ": (" + curvePoints[i].getX() + ", " + curvePoints[i].getY() + ")");
        }
        System.out.println("velo X: \n");
        for (int i = 0; i < velocityVectors.size(); i++) {
            Vector2d vector = (Vector2d) velocityVectors.get(i);
            System.out.println(vector.getX());
        }
        System.out.println("velo Y: \n");
        for (int i = 0; i < velocityVectors.size(); i++) {
            Vector2d vector = (Vector2d) velocityVectors.get(i);
            System.out.println(vector.getY());
        }
        System.out.println("velo Mag: \n");
        for (int i = 0; i < velocityVectors.size(); i++) {
            Vector2d vector = (Vector2d) velocityVectors.get(i);
            System.out.println(vector.length());
        }
        System.out.println("cent Force: \n");
        for (double force: velocityMagnitude) {
            System.out.println(force);
        }

        for (int i = 0; i < circleCenters.length-1; i++) {
            System.out.println("Point " + i + ": (" + circleCenters[i].getX() + ", " + circleCenters[i].getY() + ")");
        }
    }
}