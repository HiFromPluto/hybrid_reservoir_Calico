package bsim.d1;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

/** Minimal little-endian float64 {@code .npy}/{@code .npz} reader for D1 fixtures. */
public final class NpzReader {

    private NpzReader() { }

    public static Map<String, Array> loadNpz(Path path) throws IOException {
        Map<String, Array> out = new HashMap<>();
        try (InputStream raw = java.nio.file.Files.newInputStream(path);
             ZipInputStream zip = new ZipInputStream(raw)) {
            ZipEntry entry;
            while ((entry = zip.getNextEntry()) != null) {
                if (entry.isDirectory()) {
                    continue;
                }
                String name = entry.getName();
                if (name.endsWith(".npy")) {
                    name = name.substring(0, name.length() - 4);
                }
                out.put(name, loadNpy(readAll(zip)));
                zip.closeEntry();
            }
        }
        return out;
    }

    public static Array loadNpy(byte[] bytes) throws IOException {
        ByteArrayInputStream in = new ByteArrayInputStream(bytes);
        byte[] magic = in.readNBytes(6);
        if (magic.length < 6 || magic[0] != (byte) 0x93 || magic[1] != 'N' || magic[2] != 'U'
                || magic[3] != 'M' || magic[4] != 'P' || magic[5] != 'Y') {
            throw new IOException("not an NPY file");
        }
        int major = in.read();
        int minor = in.read();
        int headerLen;
        if (major == 1) {
            byte[] hl = in.readNBytes(2);
            headerLen = ByteBuffer.wrap(hl).order(ByteOrder.LITTLE_ENDIAN).getShort() & 0xffff;
        } else if (major == 2 || major == 3) {
            byte[] hl = in.readNBytes(4);
            headerLen = ByteBuffer.wrap(hl).order(ByteOrder.LITTLE_ENDIAN).getInt();
        } else {
            throw new IOException("unsupported NPY version " + major + "." + minor);
        }
        String header = new String(in.readNBytes(headerLen), StandardCharsets.US_ASCII);
        Matcher descr = Pattern.compile("'descr':\\s*'([^']+)'").matcher(header);
        if (!descr.find()) {
            throw new IOException("npy descr missing");
        }
        String dtype = descr.group(1);
        Matcher fortran = Pattern.compile("'fortran_order':\\s*(True|False)").matcher(header);
        if (!fortran.find()) {
            throw new IOException("npy fortran_order missing");
        }
        boolean fortranOrder = fortran.group(1).equals("True");
        Matcher shapeM = Pattern.compile("'shape':\\s*\\(([^)]*)\\)").matcher(header);
        if (!shapeM.find()) {
            throw new IOException("npy shape missing");
        }
        int[] shape = parseShape(shapeM.group(1));
        byte[] payload = in.readAllBytes();
        if (!dtype.equals("<f8") && !dtype.equals("|f8")) {
            throw new IOException("expected <f8, got " + dtype);
        }
        int n = payload.length / 8;
        double[] data = new double[n];
        ByteBuffer buf = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        for (int i = 0; i < n; i++) {
            data[i] = buf.getDouble();
        }
        if (fortranOrder && shape.length > 1) {
            throw new IOException("fortran_order arrays are not supported");
        }
        return new Array(shape, data);
    }

    private static int[] parseShape(String inner) {
        inner = inner.trim();
        if (inner.isEmpty()) {
            return new int[0];
        }
        String[] parts = inner.split(",");
        java.util.List<Integer> dims = new java.util.ArrayList<>();
        for (String p : parts) {
            p = p.trim();
            if (p.isEmpty()) {
                continue;
            }
            dims.add(Integer.parseInt(p));
        }
        int[] shape = new int[dims.size()];
        for (int i = 0; i < dims.size(); i++) {
            shape[i] = dims.get(i);
        }
        return shape;
    }

    private static byte[] readAll(InputStream in) throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        byte[] buf = new byte[8192];
        int n;
        while ((n = in.read(buf)) >= 0) {
            out.write(buf, 0, n);
        }
        return out.toByteArray();
    }

    public static final class Array {
        public final int[] shape;
        public final double[] data;

        Array(int[] shape, double[] data) {
            this.shape = shape;
            this.data = data;
        }

        public double scalar() {
            if (data.length != 1) {
                throw new IllegalStateException("not a scalar, n=" + data.length);
            }
            return data[0];
        }

        public double[] vector() {
            return data.clone();
        }

        public double[][] matrix2() {
            if (shape.length != 2) {
                throw new IllegalStateException("not a 2d array");
            }
            int r = shape[0];
            int c = shape[1];
            double[][] m = new double[r][c];
            for (int i = 0; i < r; i++) {
                System.arraycopy(data, i * c, m[i], 0, c);
            }
            return m;
        }
    }
}
