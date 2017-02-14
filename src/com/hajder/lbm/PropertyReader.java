package com.hajder.lbm;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Property file reader.
 * Reads all required data for processing LBM.
 * @author Piot Hajder
 */
final class PropertyReader {
    private static final String PROPERTY_FILE = "init_params.properties";

    static final String DEBUG = "debug";
    static final String MAX_STEP = "max_step";
    static final String EP = "ep";
    static final String SIZE_X = "size_x";
    static final String SIZE_Y = "size_y";
    static final String REYNOLDS = "Re";
    static final String TAU = "tau";

    private static Properties properties;

    /**
     * Reads property for given key as parameter.
     * Initializes property if doesn't exist.
     * @param key Property key
     * @return Property as String
     * @throws IOException If given input stream is interrupted.
     */
    static String readProperty(String key) throws IOException {
        if(properties == null) {
            InputStream fis = new FileInputStream(PROPERTY_FILE);
            properties = new Properties();
            properties.load(fis);
            fis.close();
        }
        return properties.getProperty(key);
    }
}
