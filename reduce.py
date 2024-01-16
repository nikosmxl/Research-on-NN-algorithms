import argparse
import numpy as np
from tensorflow.keras import layers, models
from sklearn.model_selection import train_test_split

def read_mnist_images(file_path):
    with open(file_path, 'rb') as f:
        # Read the magic number, number of images, number of rows, and number of columns
        magic_number = int.from_bytes(f.read(4), 'big') # Not using it anywhere.
        num_images = int.from_bytes(f.read(4), 'big')
        num_rows = int.from_bytes(f.read(4), 'big')
        num_cols = int.from_bytes(f.read(4), 'big')

        # Read image data
        image_size = num_rows * num_cols
        images = np.frombuffer(f.read(image_size * num_images), dtype=np.uint8)
        images = images.reshape((num_images, num_rows, num_cols))

    return images

def build_autoencoder(x, y, inChannel, num_filters, filter_size, num_layers, latent_dimension):
    autoencoder = models.Sequential()
    act_func = 'relu'
    last_act_func = 'sigmoid'

    # Encoder
    current_filters = num_filters

    for _ in range(num_layers - 1):
        autoencoder.add(layers.Conv2D(current_filters, kernel_size=(filter_size, filter_size), activation=act_func, padding='same'))
        autoencoder.add(layers.MaxPooling2D((2, 2), padding='same'))

        # Double the number of filters for the next layer
        current_filters *= 2

    autoencoder.add(layers.Conv2D(current_filters, kernel_size=(filter_size, filter_size), activation=act_func, padding='same'))

    encoder_output_shape = (x // 2**(num_layers - 1), y // 2**(num_layers - 1), current_filters)

    autoencoder.add(layers.Flatten())
    autoencoder.add(layers.Dense(latent_dimension, activation=act_func))

    # Decoder
    autoencoder.add(layers.Dense(np.prod(encoder_output_shape), activation=act_func))
    autoencoder.add(layers.Reshape(encoder_output_shape))

    for _ in range(num_layers - 1):
        autoencoder.add(layers.Conv2D(current_filters, kernel_size=(filter_size, filter_size), activation=act_func, padding='same'))
        autoencoder.add(layers.UpSampling2D((2, 2)))

        # Halve the number of filters for the next layer
        current_filters //= 2

    autoencoder.add(layers.Conv2D(1, kernel_size=(filter_size, filter_size), activation=last_act_func, padding='same'))

    return autoencoder

def save_encoded_data(encoder, images, output_file):
    encoded_data = encoder.predict(images)
    encoded_data_int = (encoded_data * 255).astype(np.uint8)

    # Write to a binary file with a custom header which our read functions ignore
    header = 'arbitrary words'.ljust(16, '\x00').encode('utf-8')

    with open(output_file, 'wb') as f:
        f.write(header)
        f.write(encoded_data_int.tobytes())

def main():
    parser = argparse.ArgumentParser(description='Reduce dimensionality of images using autoencoder.')
    parser.add_argument('-d', '--input-file', required=True)
    parser.add_argument('-od', '--output-dataset-file', required=True)
    parser.add_argument('-q', '--query-file', required=True)
    parser.add_argument('-oq', '--output-query-file', required=True)

    args = parser.parse_args()

    # File path to the binary file containing MNIST images
    input_file = args.input_file
    query_file = args.query_file

    # Read MNIST images into a NumPy array
    input_data = read_mnist_images(input_file)
    query_data = read_mnist_images(query_file)

    x, y = input_data.shape[1:]
    inChannel = 1

    # Rescale pixel values to be in [0,1]
    input_data = input_data.astype('float32') / 255.0
    query_data = query_data.astype('float32') / 255.0

    # Reshape the images to (28, 28, 1) for convolutional layers
    input_data = input_data.reshape((len(input_data), x, y, inChannel))
    query_data = query_data.reshape((len(query_data), x, y, inChannel))

    # Split the dataset into training and validation sets
    x_train, x_val = train_test_split(input_data, test_size=0.1, random_state=13)

    # Define hyperparameters
    num_filters = 16
    filter_size = 3
    num_layers = 3
    num_epochs = 5
    batch_size = 128
    latent_dimension = 10

    # Build the autoencoder model
    autoencoder = build_autoencoder(x, y, inChannel, num_filters, filter_size, num_layers, latent_dimension)

    # Compile the model
    autoencoder.compile(optimizer='adam', loss='mean_squared_error')

    # Train the autoencoder
    autoencoder.fit(x_train, x_train,
                    epochs=num_epochs,
                    batch_size=batch_size,
                    shuffle=True,
                    validation_data=(x_val, x_val))
    
    # Evaluate the autoencoder on the validation set
    evaluation_loss = autoencoder.evaluate(x_val, x_val)
    print(f'Validation Loss: {evaluation_loss}')

    # Extract the latent representations using the encoder part
    encoder = models.Sequential()

    for i in range(2 * num_layers - 1 + 2):  # 2 * num_layers - 1 is the convolutional/maxpooling layers + 2 for flatten and dense layer
        encoder.add(autoencoder.layers[i])

    # Save the encoded dataset to the output file
    save_encoded_data(encoder, input_data, args.output_dataset_file)

    # Save the encoded queryset to the output file
    save_encoded_data(encoder, query_data, args.output_query_file)

if __name__ == "__main__":
    main()
