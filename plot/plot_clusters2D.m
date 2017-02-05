
function plot_clusters2D(points_file, clusters_file, centroids_file)

X = load(points_file);
X(1,:) = []; % remove metadata

clusters = load(clusters_file);

centroids = load(centroids_file);
centroids(1,:) = []; % remove metadata

f = figure();
hold on

nclusters = size(centroids,1);
for c=1:nclusters
  points = X(clusters == c-1, :);
  scatter(points(:,1), points(:,2));

  scatter(centroids(c,1), centroids(c,2), 'k');%, 'LineWidth', 5, 'SizeData', 100);
end

uiwait(f);
