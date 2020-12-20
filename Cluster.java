import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.ListIterator;
import java.util.Random;
import java.util.Scanner;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

public class Cluster {
	private static final String DRAWINGS_PATH = "./full_simplified_parachute.ndjson";
	private static final int NUMBER_OF_PICTURES = 1000;
	private static final int m = 50;
	private static final int K = 7;
	private static final String OUTPUT_PATH = "./output.txt";
	private static final String MEAN1_PATH = "./mean_1.txt";
	private static final String MEAN2_PATH = "./mean_2.txt";

	public static void main(String[] args) {
		Double[] xFeatureVector = new Double[m * NUMBER_OF_PICTURES];
		Double[] yFeatureVector = new Double[m * NUMBER_OF_PICTURES];
		try {
			readFile(DRAWINGS_PATH, xFeatureVector, yFeatureVector);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (ParseException e) {
			e.printStackTrace();
		}

//		for (int i = 50; i < 100; i++)
//			System.out.print("" + xFeatureVector[i] + "," + yFeatureVector[i] + ",");
		System.out.println("done processing images");

		// generate random clusters to initialize
		Double[][] xClusters = new Double[K][m];
		Double[][] yClusters = new Double[K][m];
		for (int i = 0; i < K; i++)
			for (int j = 0; j < m; j++) {
				Random rng = new Random();
				xClusters[i][j] = (double) rng.nextInt(256);
				yClusters[i][j] = (double) rng.nextInt(256);
			}

		// as long as cluster mean changes, keep looping
		int[] kStar;
		boolean changed = false;
		do {
			// get closest clusters
			kStar = getClosestClusters(xClusters, yClusters, xFeatureVector, yFeatureVector);
			// update cluster means
			changed = updateClusterMeans(kStar, xClusters, yClusters, xFeatureVector, yFeatureVector);
		} while (changed);
		System.out.println("done updating cluster means");

		try {
			// write output.txt
			BufferedWriter bwOutput = new BufferedWriter(new FileWriter(OUTPUT_PATH));
			for (int i = 0; i < NUMBER_OF_PICTURES; i++) {
				bwOutput.write("" + kStar[i]);
				bwOutput.newLine();
			}
			bwOutput.close();

			// write mean1 and mean2
			BufferedWriter mean1Out = new BufferedWriter(new FileWriter(MEAN1_PATH));
			BufferedWriter mean2Out = new BufferedWriter(new FileWriter(MEAN2_PATH));
			for (int i = 0; i < m - 1; i++) {
				mean1Out.write("" + xClusters[0][i] + "," + yClusters[0][i] + ",");
				mean2Out.write("" + xClusters[1][i] + "," + yClusters[1][i] + ",");
			}
			mean1Out.write("" + xClusters[0][m - 1] + "," + yClusters[0][m - 1]);
			mean2Out.write("" + xClusters[1][m - 1] + "," + yClusters[1][m - 1]);
			mean1Out.close();
			mean2Out.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// return whether cluster mean changed
	public static boolean updateClusterMeans(int[] kStar, Double[][] xClusters, Double[][] yClusters,
			Double[] xFeatureVector, Double[] yFeatureVector) {
		Double xCumulative[][] = new Double[K][m];
		Double yCumulative[][] = new Double[K][m];
		boolean changed = false;
		for (int i = 0; i < K; i++)
			for (int j = 0; j < m; j++) {
				xCumulative[i][j] = 0.;
				yCumulative[i][j] = 0.;
			}

		for (int i = 0; i < K; i++) {
			int nk = 0;
			for (int j = 0; j < NUMBER_OF_PICTURES; j++)
				if (kStar[j] == i) {
					nk++;
					for (int k = 0; k < m; k++) {
						xCumulative[i][k] += xFeatureVector[m * j + k];
						yCumulative[i][k] += yFeatureVector[m * j + k];
					}
				}
			// normalize by nk and round
			for (int j = 0; j < m; j++) {
				xCumulative[i][j] = (double) Math.round(xCumulative[i][j] / nk);
				yCumulative[i][j] = (double) Math.round(yCumulative[i][j] / nk);
				if (Math.round(xClusters[i][j]) != Math.round(xCumulative[i][j])
						|| Math.round(yClusters[i][j]) != Math.round(yCumulative[i][j])) {
					changed = true;
					xClusters[i][j] = xCumulative[i][j];
					yClusters[i][j] = yCumulative[i][j];
				}
			}
		}
		return changed;
	}

	public static int[] getClosestClusters(Double[][] xClusters, Double[][] yClusters, Double[] xFeatureVector,
			Double[] yFeatureVector) {

		int[] kStar = new int[NUMBER_OF_PICTURES];

		for (int i = 0; i < NUMBER_OF_PICTURES; i++) {

			Double lowestDistance = 9999999999999.;

			for (int k = 0; k < K; k++) {
				Double d = 0.;
				for (int j = 0; j < m; j++)
					d += Math.pow(xFeatureVector[m * i + j] - xClusters[k][j], 2)
							+ Math.pow(yFeatureVector[m * i + j] - yClusters[k][j], 2.);
				d = Math.sqrt(d);
				if (d <= lowestDistance) {
					lowestDistance = d;
					kStar[i] = k;
				}
			}
		}

		return kStar;
	}

	// return number of pictures
	public static int readFile(String path, Double[] xFeatureVector, Double[] yFeatureVector)
			throws FileNotFoundException, ParseException {
		File file = new File(path);
		Scanner scanner = new Scanner(file);
		int n = 0;
		for (int i = 0; i < NUMBER_OF_PICTURES && scanner.hasNextLine(); i++) {
			n = i + 1;
			JSONArray data = (JSONArray) ((JSONObject) new JSONParser().parse(scanner.nextLine())).get("drawing");
			ArrayList<Double> xCoordinates = new ArrayList<Double>();
			ArrayList<Double> yCoordinates = new ArrayList<Double>();

			for (int j = 0; j < data.size(); j++) {

				JSONArray xVals = (JSONArray) ((JSONArray) data.get(j)).get(0);
				JSONArray yVals = (JSONArray) ((JSONArray) data.get(j)).get(1);

				ListIterator<?> iter = xVals.listIterator();
				while (iter.hasNext())
					xCoordinates.add(Double.valueOf(iter.next().toString()));

				iter = yVals.listIterator();
				while (iter.hasNext())
					yCoordinates.add(Double.valueOf(iter.next().toString()));

			}

			// upsample/downsample
			if (xCoordinates.size() > m) {
				Double[] xDownSampled = downSample(xCoordinates);
				Double[] yDownSampled = downSample(yCoordinates);
				for (int j = 0; j < m; j++) {
					xFeatureVector[j + i * m] = xDownSampled[j];
					yFeatureVector[j + i * m] = yDownSampled[j];
				}
			} else if (xCoordinates.size() < m) {
				Double[] xUpSampled = upSample(xCoordinates);
				Double[] yUpSampled = upSample(yCoordinates);
				for (int j = 0; j < m; j++) {
					xFeatureVector[j + i * m] = xUpSampled[j];
					yFeatureVector[j + i * m] = yUpSampled[j];
				}
			} else { // no upsampling or downsampling necessary?
				for (int j = 0; j < m; j++) {
					xFeatureVector[j + i * m] = xCoordinates.get(j);
					yFeatureVector[j + i * m] = yCoordinates.get(j);
				}
			}

		}

		scanner.close();
		return n;
	}

	public static Double[] downSample(ArrayList<Double> coordinates) {
		int t = coordinates.size();
		Double[] downSampled = new Double[m];
		for (int i = 0; i < m; i++)
			downSampled[i] = coordinates.get(Math.round(i * (t - 1) / (m - 1)));
		return downSampled;
	}

	public static Double[] upSample(ArrayList<Double> coordinates) {
		int t = coordinates.size();
		int mm = m;

		// check if will need to down sample
		boolean needToDownSample = false;
		if ((m - 1) % (t - 1) != 0) {
			needToDownSample = true;
			mm = (int) Math.ceil((m - 1.) / (t - 1.)) * (t - 1) + 1;
		}
		int q = (mm - 1) / (t - 1);
		Double[] upSampled = new Double[mm];

		for (int i = 0; i < mm - 1; i++) {
			upSampled[i] = (1 - ((double) i % q) / (double) q)
					* coordinates.get((int) Math.floor((double) i / (double) q));
			upSampled[i] += ((double) i % q) / (double) q
					* coordinates.get((int) Math.floor((double) i / (double) q) + 1);
		}
		upSampled[mm - 1] = coordinates.get(t - 1);

		if (needToDownSample) {
			Double[] tmp = new Double[m];
			for (int i = 0; i < m; i++)
				tmp[i] = upSampled[(int) Math.round(i * (mm - 1.) / (m - 1.))];
			upSampled = tmp;
		}

		return upSampled;
	}

}
