from backend.utils.data_generator import PolymerDataGenerator, save_dataset
from backend.utils.train_model import main as train_model

def main():
    print("Generating synthetic polymer dataset...")
    generator = PolymerDataGenerator()
    dataset = generator.generate_dataset(num_samples=10000)
    save_dataset(dataset)
    print("Dataset generated and saved!")
    
    print("\nTraining GNN model...")
    train_model()
    print("Training complete!")

if __name__ == "__main__":
    main() 