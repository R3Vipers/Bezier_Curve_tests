package org.example;

import java.awt.geom.Point2D;
import java.util.ArrayList;


public class BezierCurve2 {
    static int numPoints = 100;
    static Point2D.Double[] controlPoints = {
            new Point2D.Double(0.0, 0.0),
            new Point2D.Double(5.00, 10.00),
            new Point2D.Double(10.00, -5.00),
            new Point2D.Double(15.00, 0.0)
    };
    // Generate points on a Bézier curve
    public static Point2D.Double[] generateBezierCurve() {
        Point2D.Double[] curvePoints = new Point2D.Double[numPoints];

        for (int j = 0; j < numPoints; j++) {
            double t = (double) j / (numPoints);

            double x = (Math.pow((1 - t), 3) * controlPoints[0].getX()) +
                       (3 * Math.pow((1 - t), 2) * t * controlPoints[1].getX()) +
                       (3 * (1 - t) * Math.pow(t,2) * controlPoints[2].getX()) +
                       (Math.pow(t,3) * controlPoints[3].getX());
            double y = (Math.pow((1 - t), 3) * controlPoints[0].getY()) +
                       (3 * Math.pow((1 - t), 2) * t * controlPoints[1].getY()) +
                       (3 * (1 - t) * Math.pow(t,2) * controlPoints[2].getY()) +
                       (Math.pow(t,3) * controlPoints[3].getY());

            curvePoints[j] = new Point2D.Double(x, y);
        }

        return curvePoints;
    }
    public static ArrayList <Vector2d> generateVelocityVector(){
        ArrayList <Vector2d> velocityVectors = new ArrayList<>();
        for(int i = 0; i < numPoints; i++) {
            double t = (double) i / (numPoints - 1);
            double dxdt = firstDerivativeX(t);
            double dydt = firstDerivativeY(t);
            Vector2d vector = new Vector2d(dxdt, dydt);
            velocityVectors.add(vector);
        }

        return velocityVectors;
    }
    public static double [] generateCentripetalForceVectorMagnitude( double mass, double velocity) {
        double [] centripetalForceVectorMagnitude = new double[numPoints-1];
        for(int i = 0; i < numPoints-1; i++) {
            double t = (double) i / (numPoints);
            double dxdt = firstDerivativeX(t);
            double dydt = firstDerivativeY(t);
            double d2xdt2 = secondDerivativeX(t);
            double d2ydt2 = secondDerivativeY(t);
            double r = radius(dxdt, dydt, d2xdt2, d2ydt2);
            double force = mass*((Math.pow(velocity, 2))/r);
            centripetalForceVectorMagnitude [i] = force;
        }
        return centripetalForceVectorMagnitude;
    }
    public static Point2D.Double[] generateCircleCenter(Point2D.Double[] curvePoints) {
        Point2D.Double[] circleCenter = new Point2D.Double[curvePoints.length-1];

        for(int i = 0; i < curvePoints.length-1; i++) {
            double t = (double) i / (numPoints);
            double dxdt = firstDerivativeX(t);
            double dydt = firstDerivativeY(t);
            double d2xdt2 = secondDerivativeX(t);
            double d2ydt2 = secondDerivativeY(t);
            double r = radius(dxdt, dydt, d2xdt2, d2ydt2);
            double dxtdyt = -dxdt/dydt;
            double x = curvePoints[i].getX()- (r * Math.cos(Math.atan(dxtdyt)));
            double y = curvePoints[i].getY()+ (r * Math.sin(-Math.atan(dxtdyt)));
            circleCenter[i] = new Point2D.Double(x, y);
        }
        return  circleCenter;
    }
    public static double firstDerivativeX (double t) {
        return 3 * (((Math.pow((1-t),2)) * (controlPoints[1].getX()-controlPoints[0].getX())) +
                (2*(1-t) * t * (controlPoints[2].getX()-controlPoints[1].getX()) +
                        (Math.pow(t,2) * (controlPoints[3].getX()-controlPoints[2].getX()))));
    }
    public static double firstDerivativeY (double t) {
        return 3.0 * (((Math.pow((1-t),2)) * (controlPoints[1].getY()-controlPoints[0].getY())) +
                (2.0*(1-t) * t * (controlPoints[2].getY()-controlPoints[1].getY())+
                        (Math.pow(t,2) * (controlPoints[3].getY()-controlPoints[2].getY()))));
    }
    public static double secondDerivativeX(double t) {
        return -6 * (controlPoints[1].getX()-controlPoints[0].getX()) * (1-(t)) +
                6 * (controlPoints[2].getX()-controlPoints[1].getX()) * (1-(2*t)) +
                6 * (t) * (controlPoints[3].getX()-controlPoints[2].getX());
    }
    public static double secondDerivativeY(double t) {
        return -6 * (controlPoints[1].getY()-controlPoints[0].getY()) * (1-(t)) +
                6 * (controlPoints[2].getY()-controlPoints[1].getY()) * (1-(2*t)) +
                6 * (t) * (controlPoints[3].getY()-controlPoints[2].getY());
    }
    public static double radius (double dxdt, double dydt, double d2xdt2, double d2ydt2) {
        double mper = -dxdt/dydt;
        double dsec = (((dxdt*d2ydt2)-(dydt*d2xdt2))/Math.pow(dxdt, 2.0))/dxdt;
        return -(Math.pow(Math.pow(1 + Math.pow(mper, 2.0), 0.5)/mper,3.0)/dsec);
    }
    // Example usage
    public static void main(String[] args) {
        double mass = 11.3398;
        double velocity = 1.2;
        Point2D.Double[] curvePoints = generateBezierCurve();
        ArrayList<Vector2d> velocityVectors = generateVelocityVector();
        double [] centripetalForces = generateCentripetalForceVectorMagnitude(mass, velocity);
        Point2D.Double[] circleCenters = generateCircleCenter(curvePoints);
        // Print the generated points
        System.out.println("points: \n");
        for (int i = 0; i < numPoints; i++) {
            System.out.println("Point " + i + ": (" + curvePoints[i].getX() + ", " + curvePoints[i].getY() + ")");
        }
        System.out.println("velocity X: \n");
        for (Vector2d velocityVector : velocityVectors) {
            System.out.println(velocityVector.getX());
        }
        System.out.println("velocity Y: \n");
        for (Vector2d velocityVector : velocityVectors) {
            System.out.println(velocityVector.getY());
        }
        System.out.println("velocity Mag: \n");
        for (Vector2d velocityVector : velocityVectors) {
            System.out.println(velocityVector.length());
        }
        System.out.println("cent Force: \n");
        for (double force: centripetalForces) {
            System.out.println(force);
        }
        for (int i = 0; i < circleCenters.length-1; i++) {
            System.out.println("Point " + i + ": (" + circleCenters[i].getX() + ", " + circleCenters[i].getY() + ")");
        }
    }
}