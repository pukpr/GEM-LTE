import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import csv

def parse_sections(raw_text):
    lines = raw_text.strip().splitlines()
    data = {}
    in_section = False
    section_data = []
    section_index = 0
    current_label = None

    for i, line in enumerate(lines):
        line = line.strip()

        # Check for custom label line
        if line.startswith("[] "):
            current_label = line[3:].strip()

        # Start of section
        elif "tidal" in line.lower():
            in_section = True
            section_data = []

        # End of section
        elif "lte" in line.lower() and in_section:
            section_index += 1
            label = current_label if current_label else f"Category {section_index}"
            data[label] = section_data
            current_label = None
            in_section = False

        # Read CSV line within a section
        elif in_section and line:
            try:
                reader = csv.reader([line])
                parsed = next(reader)
                if len(parsed) >= 2:
                    label = parsed[0].strip()[:9]  # Limit to 9 characters
                    amplitude = abs(float(parsed[1].strip()))
                    section_data.append((label, amplitude))
            except Exception:
                continue
    return data

def plot_data():
    raw_text = text_input.get("1.0", tk.END)
    sections = parse_sections(raw_text)

    if not sections:
        messagebox.showerror("Error", "No valid data found.")
        return

    all_labels = sorted(set(label for pairs in sections.values() for label, _ in pairs))
    categories = list(sections.keys())
    fig, ax = plt.subplots(figsize=(10, 6))

    bar_height = 0.8 / len(categories)
    index = list(range(len(all_labels)))

    color_cycle = plt.cm.get_cmap('tab10', len(categories))

    for i, category in enumerate(categories):
        data_dict = dict(sections[category])
        heights = [data_dict.get(label, 0) for label in all_labels]
        bar_positions = [x + i * bar_height for x in index]
        ax.barh(bar_positions, heights, bar_height, label=category, color=color_cycle(i))

    ax.set_yticks([x + bar_height * (len(categories) - 1) / 2 for x in index])
    ax.set_yticklabels(all_labels)
    ax.set_xlabel("Amplitude")
    ax.set_ylabel("Tidal Factor (days)")
    ax.set_title("Tidal Factor contribution to forcing")
    ax.legend()
    ax.invert_yaxis()

    canvas.figure = fig
    canvas.draw()

# GUI Setup
root = tk.Tk()
root.title("Multi-Section Horizontal Bar Chart")

frame = ttk.Frame(root, padding=10)
frame.pack(fill="both", expand=True)

text_input = tk.Text(frame, height=15, width=80)
text_input.pack(pady=5)

button_frame = ttk.Frame(frame)
button_frame.pack(pady=5)

plot_button = ttk.Button(button_frame, text="Generate Chart", command=plot_data)
plot_button.pack(side="left", padx=5)

exit_button = ttk.Button(button_frame, text="Exit", command=root.quit)
exit_button.pack(side="left", padx=5)

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=frame)
canvas.get_tk_widget().pack(fill="both", expand=True)

root.mainloop()
